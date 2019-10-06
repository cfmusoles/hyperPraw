PROCESSES="$1"
BATCH_SIZE="$2"

mpirun -np $PROCESSES --oversubscribe hyperPraw -h resources/2cubes_sphere.mtx.hgr -i 100 -m 1200 -p rHDRF -t 0 -n test_default -P -g $BATCH_SIZE -H

mpirun -np $PROCESSES --oversubscribe hyperPraw -h resources/sat14_itox_vc1130.cnf.dual.hgr -i 100 -m 1200 -p rHDRF -t 0 -n test_default -P -g $BATCH_SIZE -H

mpirun -np $PROCESSES --oversubscribe hyperPraw -h resources/ABACUS_shell_hd.mtx.hgr -i 100 -m 1200 -p rHDRF -t 0 -n test_default -P -g $BATCH_SIZE -H

mpirun -np $PROCESSES --oversubscribe hyperPraw -h resources/sparsine.mtx.hgr -i 100 -m 1200 -p rHDRF -t 0 -n test_default -P -g $BATCH_SIZE -H

mpirun -np $PROCESSES --oversubscribe hyperPraw -h resources/pdb1HYS.mtx.hgr -i 100 -m 1200 -p rHDRF -t 0 -n test_default -P -g $BATCH_SIZE -H

mpirun -np $PROCESSES --oversubscribe hyperPraw -h resources/sat14_10pipe_q0_k.cnf.primal.hgr -i 100 -m 1200 -p rHDRF -t 0 -n test_default -P -g $BATCH_SIZE -H

mpirun -np $PROCESSES --oversubscribe hyperPraw -h resources/sat14_E02F22.cnf.hgr -i 100 -m 1200 -p rHDRF -t 0 -n test_default -P -g $BATCH_SIZE -H

mpirun -np $PROCESSES --oversubscribe hyperPraw -h resources/webbase-1M.mtx.hgr -i 100 -m 1200 -p rHDRF -t 0 -n test_default -P -g $BATCH_SIZE -H

exit

mpirun -np $PROCESSES --oversubscribe hyperPraw -h resources/atmosmodj.mtx.hgr -i 30 -m 1200 -p rHDRF -t 0 -n test_default -P

mpirun -np $PROCESSES --oversubscribe hyperPraw -h resources/kkt_power.mtx.hgr -i 30 -m 1200 -p rHDRF -t 0 -n test_default -P

mpirun -np $PROCESSES --oversubscribe hyperPraw -h resources/sat14_velev-vliw-uns-2.0-uq5.cnf.dual.hgr -i 30 -m 1200 -p rHDRF -t 0 -n test_default -P


