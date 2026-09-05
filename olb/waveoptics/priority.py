"""The process priority boost of a long fidelity-2 run.

A Windows process that has no console window runs under power throttling
(EcoQoS): the scheduler parks it on the efficiency cores and it lowers the
clock. That is the state of a run that starts over ssh or through WMI, and
it is the cause of the 11 percent load that a 16-worker pool showed
(2026-09-05, `validation/campaign_resources/`). The opt-out from that
throttling is what recovers the speed. The priority class is a smaller
lever.

Windows does not pass the opt-out to a spawned child, so a process pool
must call the boost in EVERY worker. `Campaign.run(boost=True)` does that
through its pool initializer. A thread inherits the state of its process,
so the threaded route needs the parent boost only.

The default class is Above Normal. Do NOT use High: a 16-worker pool at
High starves sshd and the VS Code server, and a Remote-SSH connection then
times out (seen 2026-09-05).

Pure ctypes, no dependency. Off Windows the boost is a no-op. It never
raises.
"""

import sys

ABOVE_NORMAL_PRIORITY_CLASS = 0x00008000
HIGH_PRIORITY_CLASS = 0x00000080
# ProcessPowerThrottling of SetProcessInformation, winbase.h.
_PROCESS_POWER_THROTTLING = 4
_PROCESS_POWER_THROTTLING_EXECUTION_SPEED = 0x1


def boost_process_priority(high=False):
    """Raise the priority class of THIS process and opt it out of EcoQoS.

    Args:
        high: False (the default) sets Above Normal. True sets High, which
              can starve the ssh server; see the module docstring.

    Returns:
        True when both calls ran and Windows accepted them. False off
        Windows or when a call failed.
    """
    if not sys.platform.startswith("win"):
        return False
    try:
        import ctypes
        from ctypes import wintypes
        k32 = ctypes.WinDLL("kernel32", use_last_error=True)
        k32.GetCurrentProcess.restype = wintypes.HANDLE
        k32.SetPriorityClass.argtypes = [wintypes.HANDLE, wintypes.DWORD]
        k32.SetPriorityClass.restype = wintypes.BOOL
        hproc = k32.GetCurrentProcess()
        cls = HIGH_PRIORITY_CLASS if high else ABOVE_NORMAL_PRIORITY_CLASS
        ok_class = bool(k32.SetPriorityClass(hproc, cls))

        class _PowerThrottlingState(ctypes.Structure):
            _fields_ = [("Version", wintypes.ULONG),
                        ("ControlMask", wintypes.ULONG),
                        ("StateMask", wintypes.ULONG)]

        k32.SetProcessInformation.argtypes = [wintypes.HANDLE, ctypes.c_int,
                                              ctypes.c_void_p, wintypes.DWORD]
        k32.SetProcessInformation.restype = wintypes.BOOL
        # ControlMask names the bit, StateMask 0 clears it: the opt-out.
        state = _PowerThrottlingState(
            1, _PROCESS_POWER_THROTTLING_EXECUTION_SPEED, 0)
        ok_throttle = bool(k32.SetProcessInformation(
            hproc, _PROCESS_POWER_THROTTLING, ctypes.byref(state),
            ctypes.sizeof(state)))
        return ok_class and ok_throttle
    except Exception:
        return False


def priority_class():
    """Give the priority class of THIS process, or None off Windows."""
    if not sys.platform.startswith("win"):
        return None
    import ctypes
    from ctypes import wintypes
    k32 = ctypes.WinDLL("kernel32", use_last_error=True)
    k32.GetCurrentProcess.restype = wintypes.HANDLE
    k32.GetPriorityClass.argtypes = [wintypes.HANDLE]
    k32.GetPriorityClass.restype = wintypes.DWORD
    return int(k32.GetPriorityClass(k32.GetCurrentProcess()))


if __name__ == "__main__":
    ok = boost_process_priority()
    cls = priority_class()
    print(f"boost ok: {ok}, priority class: {cls}")
    if sys.platform.startswith("win"):
        assert ok, "the boost failed on Windows"
        assert cls == ABOVE_NORMAL_PRIORITY_CLASS, cls
    else:
        assert ok is False and cls is None
    print("priority self-check OK")
