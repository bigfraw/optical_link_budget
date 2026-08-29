"""A small thread pool for independent, CPU-bound trials.

The wave-optics trials are independent: each one makes its own phase screens
and moves its own field. So a set of trials runs across threads with no shared
state. The numpy and scipy FFT of the split step releases the GIL while it
runs, so the threads give a real speed-up on a CPU-bound trial. See the CPython
documentation of the global interpreter lock, and the numpy note that the
pocketfft transform releases it.

The class is GENERAL. It runs any set of independent callables and it keeps the
results in the input order. It is not tied to wave optics. Load it with a list
of trial callables (`run`), or map a function over a set of items (`map`).

THREAD SAFETY IS THE CALLER'S JOB. The pool runs the callables at the same
time. Each callable must touch only its own state, or read-only shared state.
A wave-optics trial meets this: it seeds its own screens, and the split step
copies the input field before it changes it.
"""

import os
from concurrent.futures import ThreadPoolExecutor, as_completed


def _default_workers():
    """Give a sensible default worker count for a CPU-bound load.

    One worker for each core, capped at 16. The FFT is single-threaded
    (pocketfft) and it releases the GIL, but the rest of a trial (the screen
    exp, the mask, the clip) holds the GIL, and the FFT itself is
    memory-bandwidth bound past a few threads. So the thread rate saturates
    well below the core count. THE CAP IS AN olb MEASUREMENT, not a hardware
    number: the scaling study of validation/waveoptics_speed/scaling_study.py
    finds the thread rate peaks at 8 to 16 workers on every case (parallel
    efficiency about 0.35 at 16 workers and about 0.15 at 32), so 16 reaches
    near the peak rate without over-subscribing. A caller that measures a
    higher optimum on its own machine passes max_workers.
    """
    return min(16, os.cpu_count() or 1)


class Threader:
    """A thread pool that runs independent callables and keeps the order.

    Attributes:
        max_workers: the number of worker threads.
    """

    def __init__(self, max_workers=None):
        """Make a threader.

        Args:
            max_workers: the number of worker threads. None takes one for each
                core, capped at 16 (see _default_workers). A value of 1 runs the
                callables one by one in the calling thread, with no pool at all.
        """
        self.max_workers = max_workers if max_workers else _default_workers()

    def map(self, fn, items, *, progress=None):
        """Run fn on each item, and give the results in the input order.

        The pool runs the calls at the same time, but the returned list follows
        the input order, not the finishing order. The FIRST exception from any
        call propagates to the caller.

        Args:
            fn:       a callable of one argument.
            items:    an iterable of arguments.
            progress: an optional callable progress(done, total). It is called
                      once for each finished item, from the calling thread,
                      so it is safe to print in it.

        Returns:
            A list of the results, one for each item, in the input order.
        """
        items = list(items)
        total = len(items)
        results = [None] * total
        if self.max_workers == 1:
            for i, item in enumerate(items):
                results[i] = fn(item)
                if progress is not None:
                    progress(i + 1, total)
            return results

        with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
            futures = {pool.submit(fn, item): i for i, item in enumerate(items)}
            done = 0
            for future in as_completed(futures):
                results[futures[future]] = future.result()
                done += 1
                if progress is not None:
                    progress(done, total)
        return results

    def run(self, jobs, *, progress=None):
        """Run a list of zero-argument callables, and keep the order.

        This is the "load it with trials" form: give a list of thunks, one for
        each trial, and get the results back in the same order.

        Args:
            jobs:     an iterable of callables that take no argument.
            progress: an optional callable progress(done, total).

        Returns:
            A list of the results, in the input order.
        """
        return self.map(lambda job: job(), jobs, progress=progress)


if __name__ == '__main__':
    import time

    # ---- 1. the order is the input order, not the finishing order ----
    # A job that sleeps LONGER for a SMALLER index finishes in reverse, so an
    # order-preserving map must still return [0, 1, 2, 3].
    def slow_square(k):
        time.sleep(0.05 * (4 - k))
        return k * k

    th = Threader(max_workers=4)
    assert th.map(slow_square, range(4)) == [0, 1, 4, 9]

    # ---- 2. the threads really run at the same time ----
    # Four 0.2 s sleeps across four workers finish in about 0.2 s, not 0.8 s.
    # The sleep releases the GIL, exactly as the FFT does.
    def nap(_):
        time.sleep(0.2)
        return 1

    t0 = time.perf_counter()
    got = Threader(max_workers=4).map(nap, range(4))
    wall = time.perf_counter() - t0
    assert sum(got) == 4
    assert wall < 0.5, wall

    # ---- 3. max_workers=1 runs in the calling thread, still in order ----
    assert Threader(max_workers=1).map(slow_square, range(4)) == [0, 1, 4, 9]

    # ---- 4. run() takes a list of thunks ----
    jobs = [(lambda v=v: v + 100) for v in range(5)]
    assert Threader(max_workers=3).run(jobs) == [100, 101, 102, 103, 104]

    # ---- 5. the progress callback fires once for each item ----
    seen = []
    Threader(max_workers=2).map(lambda k: k, range(6),
                                progress=lambda d, t: seen.append((d, t)))
    assert [d for d, _ in seen] == [1, 2, 3, 4, 5, 6], seen
    assert all(t == 6 for _, t in seen), seen

    # ---- 6. the first exception propagates ----
    def boom(k):
        if k == 2:
            raise ValueError('boom at 2')
        return k

    try:
        Threader(max_workers=3).map(boom, range(5))
        raise AssertionError('the exception must propagate')
    except ValueError as exc:
        assert 'boom at 2' in str(exc), str(exc)

    print(f"default workers            {_default_workers():9d}")
    print(f"four 0.2 s naps, 4 workers {wall:9.3f} s")
    print("self-check passed")
