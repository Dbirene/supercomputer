Инструкции по сборке и запуску программы
Данный документ содержит инструкции по сборке и запуску программ OpenMP, MPI и OpenMP+MPI на суперкомпьютере IBM Polus.

1. Программа OpenMP (task1/cd.cpp)
   
Команды для сборки:

mkdir task1

cd task1

nano cd.cpp

Скопировать код

сохранить код ctrl+O enter

Выйти из редактирования ctrl+X

module load gcc

g++ -std=c++11 -fopenmp cd.cpp -o openmp

LSF-скрипт задания (nano openmp.lsf):

#LSBATCH: User input

#BSUB -n 4

#BSUB -W 00:15

#BSUB -o "openmp.%J.out"

#BSUB -e "openmp.%J.err"

#BSUB -R "span[hosts=1]"

OMP_NUM_THREADS=1 ./openmp

Отправка задания:

bsub < openmp.lsf

Затем измените количество потоков в файле lsf, чтобы получить результаты для разного количества потоков в одной и той же сетке. Затем измените размер сетки в файле cpp и повторите описанную выше операцию.

Тогда вы получите результаты, подобные следующим:

![7923822f9f53ac804c63e41cbd25a4d](https://github.com/user-attachments/assets/921d4e01-29e2-4635-af0c-d4f0562d5ee9)

2. Программа MPI (task2/mpi.cpp)
   
Команды для сборки:

mkdir task2

cd task2

nano mpi.cpp

Скопировать код

сохранить код ctrl+O enter

Выйти из редактирования ctrl+X

mpicxx -std=c++11 -fopenmp mpi.cpp -o mpi -lm

LSF-скрипт задания (mpi.lsf):

#LSBATCH: User input

#!/bin/bash

#BSUB -J mpi

#BSUB -o mpi_output.%J

#BSUB -e mpi_error.%J

#BSUB -n 4

#BSUB -R "span[ptile=2]"

#BSUB -q normal

#BSUB -W 00:30

#BSUB -P project_code

module load OpenMPI/4.0.0

mpirun -np 4 ./mpi

Отправка задания:

bsub < mpi.lsf

Затем измените количество процессов в файле lsf, чтобы получить результаты с разным количеством процессов в одной сетке. Затем измените размер сетки в файле cpp и повторите вышеописанное.

Тогда вы получите результаты, подобные следующим:

<img width="449" alt="9b9d460830fd5646d126d39c60e0648" src="https://github.com/user-attachments/assets/39dd95cc-f183-4992-aa8c-9d3d9c0111d9">

3. Программа OpenMP+MPI (task3/openmp_mpi.cpp)
   
Команды для сборки:

mkdir task3

cd task3

nano openmp_mpi.cpp

Скопировать код

сохранить код ctrl+O enter

Выйти из редактирования ctrl+X

LSF-скрипт задания (openmp_mpi.lsf):

#LSBATCH: User input

#!/bin/bash

#BSUB -J openmp_mpi

#BSUB -o openmp_mpi_output.%J

#BSUB -e openmp_mpi_error.%J

#BSUB -q normal

#BSUB -n 2

#BSUB -R "span[ptile=4]"

#BSUB -x

module load OpenMPI/4.0.0

module load GCC/9.3.0

export OMP_NUM_THREADS=1

mpicxx -fopenmp -std=c++11 -o openmp_mpi openmp_mpi.cpp -lm

mpirun -np 2 ./openmp_mpi

Отправка задания:

bsub < task3_openmp_mpi.lsf

Затем измените количество потоков и процессов в файле lsf, чтобы получить результаты с разным количеством потоков и процессов в одной сетке. Затем измените размер сетки в файле cpp и повторите описанное выше.

Тогда вы получите результаты, подобные следующим:

<img width="416" alt="8cc0808ed4d0c832f280fe79234d290" src="https://github.com/user-attachments/assets/ff51dd3a-f58b-48c7-bb49-00b341421363">

Примечания

Измените параметры #BSUB -n и #BSUB -q в зависимости от количества необходимых узлов и типа очереди.

Вывод программы (*.out и *.err) будет сохранён в соответствующих файлах. Проверьте их содержимое при необходимости.




