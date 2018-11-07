#include <mpi.h> 
#include "stdio.h" 
#include "stdlib.h" 
int main(int argc, char** argv){ 
  double Tmax; 
  int n, m1, m2; 
  double ht, hx, hy; 
  double aT, lambda, CT, p; 
  int p_rank, t_size; 
  double *T;   // для корневого (root) процесса
  double *tT;  // для процессов - "потомков"
  int mbuf[3]; // хранит размерность массива для передачи процессам - "потомкам"
  /* Инициализация MPI */ 
  MPI_Init (&argc, &argv); 
  MPI_Status p_status;   // число процессов  
  MPI_Comm_size (MPI_COMM_WORLD, &t_size);  // ранк текущего процесса
  MPI_Comm_rank (MPI_COMM_WORLD, &p_rank); 
  printf ("Total size is %d processes, current pid-rank is %d.\n", t_size, p_rank); 
  if(!p_rank){ 
        /* Инициализация данных в корневом процессе */ 
        if(argc<2){
            m1=10; m2=10; Tmax=10.0; ht=0.1; hx=0.1; hy=0.1;
        } 
        else if(argc==2){
            m1=atoi(argv[1]); m2=m1; Tmax=atoi(argv[1]); 
            ht=0.1; hx=0.1; hy=0.1;
        } 
        else if(argc==3){
            m1=atoi(argv[1]); m2=m1; Tmax=atoi(argv[2]); 
            ht=0.1; hx=0.1; hy=0.1;
        } 
        else if(argc==4){
            m1=atoi(argv[1]); m2=atoi(argv[2]); 
            Tmax=atoi(argv[3]); ht=0.1; hx=0.1; hy=0.1;
        } 
        else if(argc==5){
            m1=atoi(argv[1]); m2=atoi(argv[2]); 
            Tmax=atoi(argv[3]); ht=atoi(argv[5]); hx=atoi(argv[5]); 
            hy=atoi(argv[5]);
        } 
        else if(argc==6){
            m1=atoi(argv[1]); m2=atoi(argv[2]); 
            Tmax=atoi(argv[3]); ht=atoi(argv[6]); hx=atoi(argv[5]); 
            hy=atoi(argv[5]);
        } 
        else if(argc>=7){
            m1=atoi(argv[1]); m2=atoi(argv[2]); 
            Tmax=atoi(argv[3]); ht=atoi(argv[7]); hx=atoi(argv[5]); 
            hy=atoi(argv[6]);} 
            ht=0.1; hx=0.1; hy=0.1; n=(int)Tmax/ht; n=Tmax/ht; 
            m1/=hx; m2/=hy; 
            lambda = 210.0; 
            CT = 880.0; 
            p = 2700.0; 
            aT=lambda/(CT*p); 
            T = (double *) malloc (sizeof(double)*n*m1*m2); 
            for (int k=0; k<n; k++) 
              for (int i=0; i<m1; i++) 
                for (int j=0; j<m2; j++){ 
                    if(k==0){//if k 
                        T[k*(m1-1)*(m2-1) + i*(m2-1) + j]=1.0; 
                        /* Граничные условия I рода */ 
                        if( (i==0) || (j==0) || (i==m1-1) || (j==m2-1) ){//left || bottom || right || top 
                            T[k*(m  1-1)*(m2-1) + i*(m  2-1) + j]=200.0; 
                        } 
                        /* Граничные условия II рода */ 
                        if( (i>0) && (j>0) && (i<=m1/2) && (j<=m2/2) ){ 
                            T[k*(m1-1)*(m2-1) + i*(m2-1) + j]=200.0*ht - T[k*(m1-1)*(m2-1) + (i-1)*(m2-1) + j]; 
                            T[k*(m1- 1)*(m2- 1) + i*(m2- 1) + j]=200.0*ht  -  T[k*(m1- 1)*(m2- 1) + i*(m2- 1) + j- 1];
                        }    
                    }//end if k 
                } 
            fprintf(stdout, "Rank %d, Root Array %d %d %d\n", p_rank, n, m1, m2); 
            mbuf[0]=n; mbuf[1]=m1; mbuf[2]=m2; 
  }   
  MPI_Bcast((void *)mbuf, 3, MPI_INT, 0, MPI_COMM_WORLD); 
  n=mbuf[0]; m1=mbuf[1]; m2=mbuf[2]; 
  fprintf(stdout, "Rank %d, Local array %d %d %d\n", p_rank , n, m1/t_size, m2); 
  tT=(double *) malloc (sizeof(double)*n*m1*m2/t_size); 
    for (int k=0; k<n; k++) 
      for (int i=0; i<m1/t_size; i++) 
        for (int j=0; j<m2; j++){ 
            tT[k*m1/t_size*m2 + i*m2+j]=0; 
    } 
  MPI_Scatter((void *)T, n*m1*m2/t_size, MPI_DOUBLE, 
      (void *)tT, n*m1*m2/t_size, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
  /* Основной вычислительный цикл */ 
  for (int k=0; k<n; k++) 
    for (int i=0; i<m1/t_size; i++) 
        for (int j=0; j<m2; j++){ 
            if( (i==0) || (j==0) || (i==m1/t_size-1) || (j==m2-1) ){ 
                    tT[(k+1)*(m1/t_size- 1)*(m2- 1) + i*(m2- 1) + j] =  aT*(  - 2*tT[k*(m1/t_size- 1)*(m2- 1) + i*(m2- 1) + j] )*ht/(hx*hx) +  tT[k*(m1/t_size- 1)*(m2- 1)  + i*(m2- 1) + j] + aT*( -  2*tT[k*(m1/t_size- 1)*(m2- 1) + i*(m2- 1) + j] )*ht/(hy*hy);
                
                    if(i!=0) 
                        tT[(k+1)*(m1/t_size- 1)*(m2- 1) + i*(m2- 1)  + j] += aT*tT[k*(m1/t_size- 1)*(m2- 1) + (i- 1)*(m2- 1) + j] *ht/(hx*hx);
                    if(j!=0) 
                        tT[(k+1)*(m1/t_size- 1)*(m2- 1) + i*(m2- 1)  + j] += aT*tT[k*(m1/t_size- 1)*(m2- 1) + i*(m2- 1) + j- 1] *ht/(hx*hx);
                    if(i!=(m1/t_size- 1)) 
                        tT[(k+1)*(m1/t_size- 1)*(m2- 1) + i*(m2- 1) + j] += aT*tT[k*(m1/t_size- 1)*(m2- 1) + (i+1)*(m2- 1) + j]  *ht/(hx*hx);
                    if(j!=(m2- 1)) 
                        tT[(k+1)*(m1/t_size- 1)*(m2- 1) +  i*(m2- 1) + j] += aT*tT[k*(m1/t_size- 1)*(m2- 1) + i*(m2- 1) + j+1]  *ht/(hx*hx);
            } else if( (i>0) && (j>0) &&  (i<m1/t_size- 1) && (j<m2- 1)  ){  
                    tT[(k+1)*(m1/t_size- 1)*(m2- 1) + i*(m2- 1) + j] =  aT*( tT[k*(m1/t_size- 1)*(m2- 1) + (i+1)*(m2- 1) + j] - 2*tT[k*(m1/t_size- 1)*(m2- 1) + i*(m2- 1) + j] + tT[k*(m1/t_size- 1)*(m2- 1) + (i- 1)*(m2- 1) + j] )*ht/(hx*hx) + tT[k*(m1/t_size- 1)*(m2- 1)  +  i*(m2- 1) + j] + aT*( tT[k*(m1/t_size- 1)*(m2- 1) + i*(m2- 1) + j+1] - 2*tT[k*(m1/t_size- 1)*(m2- 1) + i*(m2- 1) + j] + tT[k*(m1/t_size- 1)*(m2- 1) + i*(m2- 1) + j- 1] )*ht/(hy*hy);
            } 
            /* Обмен данными между процессами */ 
            if(t_size>1){// if 
                if( p_rank>0 && p_rank<t_size-1 && t_size>1) 
                    MPI_Sendrecv((void *)(tT+(k+1)*m1*m2), 
                        m1/t_size*m2, MPI_DOUBLE, p_rank+1, 0, (void *)(tT+k*m1*m2), 
                        m1/t_size*m2, MPI_DOUBLE, p_rank-1, 0, MPI_COMM_WORLD, &p_status); 
                else if(!p_rank && t_size>1)  { 
                        MPI_Send((void *)(tT+(k+1)*m1*m2), m1/t_size*m2, 
                        MPI_DOUBLE, p_rank+1, 0, MPI_COMM_WORLD); 
                } 
                else if( (p_rank == t_size-1) && t_size>1) {  
                            MPI_Recv((void *)(tT+k*m1*m2), m1/t_size*m2, 
                            MPI_DOUBLE, p_rank-1, 0, MPI_COMM_WORLD, &p_status);  
                } 
            }// end if 
        } 
  MPI_Barrier( MPI_COMM_WORLD ); 
  MPI_Gather((void *)tT, n*m1*m2/t_size, MPI_DOUBLE, 
    (void *)T, n*m1*m2/t_size, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
  MPI_Finalize(); 
  if(!p_rank){ 
            /* * * GNUPLOT * * */ 
            float x[m1][m2]; 
            float y[m1][m2]; 
            int k = 0; 
            for (int i=0; i<m1; i++){//for i 
                for (int j=0; j<m2; j++){//for j 
            x[i][j]=(i)*hx; 
            y[i][j]=(j)*hy; 
        } //end; for j 
    } //end; for i 
FILE *f_output, *gp_output; 
gp_output=fopen ("thermal_ani.plot","w"); 
if( gp_output == NULL ) printf("Cannot open the output file.\n"); 
f_output=fopen ("thermal.log","w"); 
if( f_output == NULL ) printf("Cannot open the output file.\n"); 
else{ 
    fprintf(gp_output, "set terminal X11\n"); //1: для вывода в окне X11 
    fprintf(gp_output, "#!gnuplot -persist thermal_ani.plot\n"); //2 
    fprintf(gp_output, "set term gif animate delay 100\n");//2: для печати в .gif-файл
    fprintf(gp_output, "set view map\n"); //2 
    fprintf(gp_output, "set output \"thermal_animate.gif\"\n");//2: для печати в .gif-файл
    fprintf(gp_output, "set xrange [%d:*]\n", 0); 
    fprintf(gp_output, "set yrange [%d:*]\n", 0); 
    fprintf(gp_output, "set xlabel 'x'; set ylabel 'y'; set cblabel 'T(x,y)'\n"); 
    fprintf(gp_output, " do for [i=0:%d] {\n", n-1); 
    fprintf(gp_output, "     splot 'thermal.log' index i u 1:2:3 w pm3d\n"); 
    fprintf(gp_output, "     pause 0.1\n"); 
    fprintf(gp_output, " }\n"); 
    
    for (int k=0; k<n; k++){//for k 
        for (int i=0; i<m1; i++){//for i 
            for (int j=0; j<m2; j++){//for j 
                fprintf(f_output,"%6f %6f %6f\n", x[i][j], y[i][j], T[k*(m1-1)*(m2-1) + i*(m2-1)+j]); 
                } //end; for j 
                fprintf(f_output,"\n"); 
            } //end; for i 
            fprintf(f_output,"\n"); 
        } //end; for k 
    fprintf(f_output,"\n"); 
    }
    fclose(f_output); 
    fclose(gp_output); 
    T = (double *) calloc (n*m1*m2, sizeof(double)); 
  }   
} 
