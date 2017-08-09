#2D 1 process
aprun -n 1 ./poisson meshes/square_level_1.bin  > log_convergence_2D_numprocs_1_lvl_1
aprun -n 1 ./poisson meshes/square_level_2.bin  > log_convergence_2D_numprocs_1_lvl_2
aprun -n 1 ./poisson meshes/square_level_3.bin  > log_convergence_2D_numprocs_1_lvl_3
aprun -n 1 ./poisson meshes/square_level_4.bin  > log_convergence_2D_numprocs_1_lvl_4
aprun -n 1 ./poisson meshes/square_level_5.bin  > log_convergence_2D_numprocs_1_lvl_5

#2D 5 processes
aprun -n 5 ./poisson meshes/square_level_1.bin  > log_convergence_2D_numprocs_5_lvl_1
aprun -n 5 ./poisson meshes/square_level_2.bin  > log_convergence_2D_numprocs_5_lvl_2
aprun -n 5 ./poisson meshes/square_level_3.bin  > log_convergence_2D_numprocs_5_lvl_3
aprun -n 5 ./poisson meshes/square_level_4.bin  > log_convergence_2D_numprocs_5_lvl_4
aprun -n 5 ./poisson meshes/square_level_5.bin  > log_convergence_2D_numprocs_5_lvl_5

#3D different num processes
aprun -n 15 ./poisson meshes/cube_level_1.bin  > log_convergence_3D_numprocs_15_lvl_1
aprun -n 20 ./poisson meshes/cube_level_2.bin  > log_convergence_3D_numprocs_20_lvl_2
aprun -n 30 ./poisson meshes/cube_level_3.bin  > log_convergence_3D_numprocs_30_lvl_3
aprun -n 32 ./poisson meshes/cube_level_4.bin  > log_convergence_3D_numprocs_32_lvl_4
aprun -n 32 ./poisson meshes/cube_level_5.bin  > log_convergence_3D_numprocs_32_lvl_5
