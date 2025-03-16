all runs with randomly initialised params
the timed call is the loss function
    simulate!(sys_dirty2, simulator, 10; n_threads=1, run_loggers=false)
only timing the simulate call (and perhaps ADAM execution)


clean - loss_copy_and_buf
1.138s
129.14 MiB

2 params - Zygote.gradient(loss, bond_k_alt, bond_r0_alt, sys_dirty,sys_clean_energy, simulator, neighbors)
sampling for 60s
6.523s
678.28 MiB
max GC 0.65%

2 params - Zygote.gradient(loss, bond_k_alt, bond_r0_alt, sys_dirty,sys_clean_energy, simulator, neighbors) and ADAM updates
sampling for 60s
6.626s
678.30MiB
max GC 0.77%


22062 params  - Zygote.jacobian(loss_all_bonds, params, sys_dirty,sys_clean_energy, simulator, neighbors)
sampling for 600 seconds
10.123 s
4.29GiB
max GC 17.80%


22062 params  - Zygote.jacobian(loss_all_bonds, params, sys_dirty,sys_clean_energy, simulator, neighbors) and ADAM updates
sampling only once
10.141 s
4.29GiB
max GC 18.06%



