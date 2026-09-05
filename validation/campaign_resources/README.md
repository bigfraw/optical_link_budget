# The 4000-trial fidelity-2 downlink Campaign with a resource monitor

`campaign_resources.py` runs one fidelity-2 space downlink through the
production `Campaign` store and it records the CPU and the RAM of the
machine while the trials run. Its purpose is to VERIFY that the process pool
does the parallel work.

The case: a 700 mm ground aperture with a 30 percent central obscuration, an
SMF detector, 1550 nm, a 500 km orbit at 30 deg, the standard preset, the
fixed outer scale L0 = 25 m, seed 20260905. The default is 4000 trials in
blocks of 50 (80 blocks) on 16 workers.

## Run it on the desktop

```
python -m validation.campaign_resources.campaign_resources
python -m validation.campaign_resources.campaign_resources --workers 24
python -m validation.campaign_resources.campaign_resources --threads
```

`--threads` gives the other level of parallelism (serial blocks, each one
threaded inside) for a comparison. `--smoke` is a small local check.

## Read the result

- `resources_<tag>.csv`: one row every `--sample-s` seconds:
  `t_s, cpu_pct, cores_busy, ram_used_mb, ram_avail_mb, n_python, python_rss_mb`.
- `resources_<tag>.json`: the mean and the peak of each column, the wall
  time and the trials per second.
- `resources_<tag>.png`: the busy cores, the RAM and the python.exe count
  against time, with the `workers` line drawn.

The pool is correct when `n_python` reads `workers + 1` while the blocks run
and `cores_busy` sits near `workers`. Keep `block_size <= n_trials / workers`,
or the pool has fewer blocks than processes (docs/api-waveoptics.md Section 9g).

The campaign directory (`campaign_<preset>_el<deg>/`) holds the blocks and
it is resumable: a second call with the same settings computes the missing
blocks only.

## Make the spawned processes really run over ssh

The desktop `bigfraw` runs a windowless ssh process throttled by default, and
a process pool adds its own traps. Do these steps, in this order. Each one
was verified on 2026-09-05 with this script.

1. **Chain remote commands with `;`, never `&&`.** The login shell of
   `ssh desktop` is Windows PowerShell 5.1, where `&&` is a parse error and
   nothing runs. Inside a `cmd /c "..."` string `&&` is correct, because
   that string runs in cmd.
2. **Ship the code from a PowerShell prompt, not from a bash shell.** Under
   bash the push script resolves `tar` to the Git tar and the archive step
   fails.
3. **Launch the run through WMI, so it outlives the ssh session.** A
   process started with `Start-Process` over ssh dies when the session
   closes. Create it with WMI and redirect its output to a log:

   ```
   $cmd = 'cmd /c "cd /d D:\repos\optical_link_budget && C:\Users\alexf\anaconda3\envs\olb\python.exe -u -m validation.campaign_resources.campaign_resources --workers 8 > validation\campaign_resources\run_w8.log 2>&1"'
   Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{CommandLine=$cmd}
   ```

   The WMI host owns the process, so it keeps running after the ssh session
   ends.
4. **Boost the parent process.** `main` calls `boost_process_priority()`
   first. It sets the Above Normal priority class and it opts the process
   out of power throttling (EcoQoS). The opt-out is what stops the
   throttling. Do NOT use the High class: a 16-worker pool at High starves
   sshd and the VS Code server, and a Remote-SSH connection then times out
   (2026-09-05).
5. **Boost every spawned worker too.** Windows does not pass that opt-out
   to a child process. With only the parent boosted, a 16-worker pool showed
   17 busy threads but about 11 percent in Task Manager: the workers ran
   throttled. The script replaces `campaign._init_worker` with a wrapper
   that calls the boost in each worker before the Campaign initializer.
   Threads inherit the parent state, so `--threads` needs step 4 only.
6. **Give the pool enough blocks.** The effective process count is
   `min(workers, ceil(n_trials / block_size))`. Keep
   `block_size <= n_trials / workers`.
7. **Judge the load with the right counter.** Task Manager shows
   `% Processor Utility`, which is frequency weighted. Read it next to
   `% Processor Time`:

   ```
   Get-Counter '\Processor Information(_Total)\% Processor Utility','\Processor Information(_Total)\% Processor Time'
   ```

   Utility above time means the cores boost. Utility below time means they
   throttle. The `cpu_mhz` column is not a throttling probe on the
   i9-14900HX: it reads a constant 1900. Confirm the process count with
   `n_python` in the CSV: it must read `workers + 1`.
8. **Kill the orphans before a relaunch.** A stopped ssh side can leave the
   workers alive. Stop the python processes whose command line contains
   `multiprocessing` or the script name, then check the count reads zero.
   The campaign resumes from the blocks on disk, so a kill loses at most
   the blocks in flight.
9. **Pick the worker count for the box.** 32 workers pegs all 32 threads
   (processor time 100 percent, utility 157 percent) and starves ssh
   itself. 8 workers is the chosen default on bigfraw.

Measured on 2026-09-05 (512 px grid, 9 screens, blocks of 50):

| workers | worker boost | trials per s | busy threads of 32 |
|---------|--------------|--------------|--------------------|
| 16      | no           | 6.7          | 17 (throttled)     |
| 8       | yes          | 7.8          | 12.7 mean, 18.8 peak |
| 16      | yes          | 7.2          | 25.5 mean, 31.2 peak |
| 12      | yes          | 8.5          | 17.5 mean, 24.4 peak |

The 12-worker run (stopped by hand at 60 of 80 blocks): 3000 trials in
355 s (0.118 s/trial), the fastest of the three. So on this 512 px grid
the plateau sits near 12 workers: 8 leaves throughput on the table, 16
adds memory traffic and slower efficiency-core threads for no gain.

The boosted 16-worker run: 4000 new trials in 559 s (0.140 s/trial), 17
python.exe while the blocks ran, 3.1 GB peak python working set, processor
time 81 percent and utility 123 percent. It is NOT faster than 8 workers on
this 512 px grid: the run is memory-bandwidth bound past about 8 processes
(see docs/api-waveoptics.md Section 9g), so 16 workers burn twice the
threads for the same rate. A 1024 px grid is where more workers pay.

The 8-worker run: 3200 new trials in 412 s (0.129 s/trial), 10 python.exe
at peak (8 workers, the parent and one tasklist child), 1.7 GB peak python
working set, 250 MB of blocks on disk for 4000 trials. The JSON of that run
reports `trials_per_s` as 9.7 because it divided ALL 4000 stored trials by
the wall time; the script now counts the new trials only.
