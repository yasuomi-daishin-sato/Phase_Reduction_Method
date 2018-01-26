#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define M 7000000

#define eps 1.0e-5            

double v[ 4 ], v1[ 4 ], q[ 6 ], q1[ 6 ];

double mu;

double z;
double z_sn;

//SN -> AH
//double a = 1.00;
//double c = 0.60;
//double d = 1.10;
//double b ;
//double c2 = 2.5;

//HC -> SN
//double a = 1.00;
//double c = 0.75;
//double d;
//double b = 2.10;
//double c2 = 2.8;

//HC -> AH
double a = 1.00;
double c = 1.00;
double d;
double b = 2.30;
double c2 = 3.2;


double theta = 0.0;
double H, T;
double c1[ M ][ 6 + 2 ];

double V, W, P1, P2, Q1, Q2;
double P1old, P2old, Q1old, Q2old;
double CC, CCnew, CCold;
double Tnew, Told, T1;

double Vfp, Wfp;
double frequency, frequency_old, frequency1;

double Vmax, Vmin;
double Va, ka;
double n_inf, b_inf;

FILE *fro, *fro1; 

//RHR Equation For Newton
double Wi1( double y ) { 
	double e = 0.0;
	
	e = ( a * y * y + b * y + c ) / d;
	
	return e;
}

double F0( double y ) {
	double e = 0.0;
	
	e = y - y * y * y / 3.0 - Wi1( y ) + z;
	
	return e;
}

//RHR 一階微分
double DWi1( double y ) {
	double e = 0.0;
	
	e = 2.0 * a * y / d + b / d;
	
	return e;
}

double DF0( double y ) {
	double e = 0.0;
	
	e = 1.0 - y * y - DWi1( y ) ;
	
	return e;
}

//Newton Method
void newton( double a1 ) {

	double newa1;

	for(;;) {
	
		newa1 = a1 - F0( a1 )/DF0( a1 );
		
		if( fabs( newa1 - a1 ) < eps ) {
			Vfp = newa1;
			//printf("%f\n", Vfp);
			break;
		}
		
		a1 = newa1;
	}
}


//For Runge-Kutta
//微分
double Df11(int j, double t, double *x) {
  double e = 0.0;
  
  e = c2  - c2 * x[ j ] * x[ j ] ;
  
  return e;
}

double Df12(int j, double t, double *x) {
  double e = 0.0;

  e = - c2;

  return e;
}

double Df21(int j, double t, double *x) {
  
	double e = 0.0;
	
	e = mu / c2 * (2.0 * a * x[ j ] + b);
	
	return e;
}

double Df22(int j, double t, double *x) {
  double e = 0.0;
  
  e = - mu * d / c2;
  
  return e;
}

double F(int j, double t, double *x) {
  double e = 0.0;
  if( j == 0 ) { 
	  e = c2 * ( x[ j ] - x[ j ] * x[ j ] * x[ j ] / 3.0 - x[ j + 1 ] + z );
	  
  }
  else if( j == 1 ) {
	  e = mu / c2 * ( a * x[ j - 1 ] * x[ j - 1 ] + b * x[ j - 1 ] + c - d * x[ j ] );
  }
  else if( j == 2 ) {
	  e = Df11( j - 2, t, (double *)x ) * x[ j ] + Df12(j - 2, t, (double *)x) * x[ j + 1 ];
  }
  else if( j == 3 ) {
      e = Df21( j - 3, t, (double *)x ) * x[ j - 1 ] + Df22( j - 3, t, (double *)x ) * x[ j ];
  }
  else {
	  e = 0.0;
  }
  return e;
}

double F1(int j, double t, double *x) {
	
	double e = 0.0;
	if( j == 0 ) { 
		e = - c2 * ( x[ j ] - x[ j ] * x[ j ] * x[ j ] / 3.0 - x[ j + 1 ] + z ) ;
		
	}
	else if( j == 1 ) {
		e = - mu / c2 * ( a * x[ j - 1 ] * x[ j - 1 ] + b * x[ j - 1 ] + c - d * x[ j ] );
	}
	else if( j == 2 ) {
		e = - Df11( j - 2, t, (double *)x ) * x[ j ] - Df12( j - 2, t,(double *)x ) * x[ j + 1 ];
	}
	else if( j == 3 ) {
		e = - Df21( j - 3, t, (double *)x ) * x[ j - 1 ] - Df22( j - 3, t, (double *)x ) * x[ j ];
	}
	else if( j == 4 ) {
		e = Df11( j - 4, t, (double *)x ) * x[ j ] + Df21( j - 4, t, (double *)x ) * x[ j + 1 ];
	}
	else if( j == 5 ) {
		e = Df12( j - 5, t, (double *)x ) * x[ j - 1 ] + Df22( j - 5, t, (double *)x ) * x[ j ];
	}
	else {
		e = 0.0;
    }
	return e;
}


void printout( double t ) {
	
	printf("%f ",t);
	//fprintf(fro, "%f ", t);
	//printf("%f ",phi);
	//for(j=1;j<=8*N;j++)
	for(int j = 0 ; j < 4 ; j++) {
		printf( "%f ", v[ j ] );
	//	fprintf( fro, "%f ", v[ j ] );
	}
	printf( "\n" );
	//fprintf(fro, "\n");

}

void printout1( int ii ) {

	printf( "%d ", ii );
	printf( "%f ", c1[ ii ][ 0 ] );
	
	for(int j = 1 ; j < 5; j ++ ) {
		printf( "%f ", c1[ ii ][ j ] );
	}
	printf( "\n" );

}

void printout2( int ii ) {
	
	fprintf( fro, "%f ", c1[ ii ][ 0 ] );
	printf( "%f ", c1[ ii ][ 0 ] );
	//
	for(int j = 1 ; j < 6 + 2 ; j ++) {
		fprintf(fro, "%f ", c1[ ii ][ j ] );
		printf( "%f ", c1[ ii ][ j ] );
	}
	fprintf( fro, "\n" );
    printf( "\n" );

	//fprintf( fro4, "%f %f %f\n", c[ ii ][ 0 ] / T1, c[ ii ][ 9 ], mu );

}

void runge( double t ) {
 
  double wv1[ 4 ], wv2[ 4 ], wv3[ 4 ], k1[ 4 ], k2[ 4 ], k3[ 4 ], k4[ 4 ];
  
  for(int j = 0 ; j < 4 ; j++ ) {
      k1[ j ] = H * F( j, t, (double *)v );
	  wv1[ j ] = v[ j ] + k1[ j ] / 2.0;
  }
  for(int j = 0 ; j < 4 ; j++ ) {
	  k2[ j ] = H * F( j, t + H / 2.0, (double *)wv1 );
	  wv2[ j ] = v[ j ] + k2[ j ] / 2.0;
  }
  for(int j = 0 ; j < 4 ; j++ ) {
      k3[ j ] = H * F( j, t + H / 2.0, (double *)wv2 );
      wv3[ j ] = v[ j ] + k3[ j ];
  }
  for(int j = 0 ; j < 4 ; j++ ) {
	  k4[ j ] = H * F( j, t + H, (double *)wv3 );
  }
  for(int j = 0 ; j < 4 ; j++) {
      v1[ j ] = v[ j ] + (k1[ j ] + 2.0 * k2[ j ] + 2.0 * k3[ j ] + k4[ j ]) / 6.0;
  }
}

void runge1( double t ) {
	
	double wv1[ 6 ], wv2[ 6 ], wv3[ 6 ], k1[ 6 ], k2[ 6 ], k3[ 6 ], k4[ 6 ];
	
	for( int j = 0 ; j < 6 ; j++ ) {
		k1[ j ] = H * F1( j, t, (double *)q );
		wv1[ j ] = q[ j ] + k1[ j ] / 2.0;
	}
	for( int j = 0 ; j < 6 ; j++ ) {
		k2[ j ] = H * F1( j, t + H / 2.0, (double *)wv1 );
		wv2[ j ] = q[ j ] + k2[ j ] / 2.0;
	}
	for( int j = 0 ; j < 6 ; j++ ) {
		k3[ j ] = H * F1( j, t + H / 2.0, (double *)wv2 );
		wv3[ j ] = q[ j ] + k3[ j ];
	}
	for( int j = 0 ; j < 6 ; j++ ) {
		k4[ j ] = H * F1( j, t + H, (double *)wv3 );
    }
	for( int j = 0 ; j < 6 ; j++ ) {
		q1[ j ] = q[ j ] + ( k1[ j ] + 2.0 * k2[ j ] + 2.0 * k3[ j ] + k4[ j ] ) / 6.0;
	}
}

int main(int i2, int j, double t, int count, int count1, int count2, int j1, int count3) {
	
	//fro = fopen("c:/gnuplot/binary/RHRClass1Class1PRCd2.1z0.073mu1.03.dat","w");
    //fro = fopen("c:/gnuplot/binary/RHRClass1Class1dfdz_dmu1.00z0.288.dat","w");
    //fro = fopen("c:/gnuplot/binary/RHRPRCCal_dmu0.50_2.dat","w");
    fro = fopen("c:/gnuplot/binary/RHRHC_AHPRCmu1.00d1.00.dat","w");
    //fro = fopen("c:/gnuplot/binary/RHRClass1Class1fmuI0.40.dat","w");
    //fro = fopen("c:/gnuplot/binary/RHRClass1Class2AvePRC_zmu0.10_2.dat","w");
    //fro = fopen("c:/gnuplot/binary/RHR_Bif_z_d_mu.dat","w");
    //fro = fopen("c:/gnuplot/binary/RHRdfdz_dmu1.00.dat","w");
    
	//d = 1.8;
	
	//刻み幅，可能処理時間，スパイクカウント初期．
    H = 0.01;
    T = 20000.0; //100000000.0;


	mu = 1.00; 

	frequency_old = 0.0;
	frequency1 = 0.0;

	//z = 0.288; //0.1685; //0.288; //0.073; //0.347;
    d = 0.5;
	Va = 65.0;
	ka = -25.0;

	//for( d = 0.5 ; d < 0.601 ; d = d + 0.001 ) {

    double V0 = 4.0;
    //double delta_z = 0.00001;
    //double delta_z1 = 0.00001;
    double z_0 = 0.25; //0.44;
	//z = 0.25; //0.073;
	for ( z = z_0 ; z > 0.0 ; z = z - 0.00005) {

	//for( z = 0.1612 ; z > 0.0 ; z = z - delta_z ) {
	//for( z = 1.00 ; z > 0.0 ; z = z - delta_z ) {
	//for( mu = 0.10 ; mu < 2.51 ; mu = mu + 0.005 ) {
	//for( mu = 2.5 ; mu > 0.5 ; mu = mu - 0.005 ) {
	//for( mu = 2.5 ; mu > 0.0 ; mu = mu - 0.005 ) {
	//for( d = 2.1 ; d > 1.8 ; d = d - 0.001 ) {
	//for( d = 1.80 ; d < 2.105 ; d = d + 0.001 ) {

	frequency = 0.0;
	  
	
	//Initial Condition Setting Up
	newton( V0 );
	Wfp = Wi1( Vfp );
	
	//printf( "%f %f\n", Vfp, Wfp);
	
	V0 = Vfp + 0.05;

	count = 0;
    
	//for( j = 0 ; j < 4 ; j++) {
	//	v[ j ] = 0.5; //(double)rand()/((double)(RAND_MAX)+1.0);
	//}
		

	v[ 0 ] = Vfp + 0.1; // -2.903123+0.1; //
	v[ 1 ] = Wfp ;  //5.470835;
	v[ 2 ] = F(0, t, (double *)v);
	v[ 3 ] = F(1, t, (double *)v);
	    
    //初期時刻
    t = 0.0;
	//printout( t );
      
    //ルンゲック処理
	for( t = H ; t < T + H ; t = t + H ) {
		runge(t);
		//if(v1[ 0 ] > Vmax ) Vmax = v1[ 0 ];
		//if(v1[ 0 ] < Vmin ) Vmin = v1[ 0 ];

		if( v[ 0 ] < theta && v1[ 0 ] > theta ) {

			count++;
			if( count == 10 ) {
			
				Told = t;
			
			}
			else if( count == 11 ) {

				Tnew = t;

			}
		  //printf("%f %d\n",t,count);
		}

		if(count == 0 && t > 500.0 && fabs(v1[ 0 ] - v[ 0 ]) < 0.00000001 ) {
			frequency = 0.0;
			break;
		}

		for(int j = 0 ; j < 4 ; j ++) {
			
			v[ j ] = v1[ j ];
		
		}
		//if(t > 700.0) printout( t );
	
		if( count == 11 ) {
			T1 = Tnew - Told;
			frequency = 1.0 / (Tnew - Told);
			
			V = v1[ 0 ];
		    W = v1[ 1 ];
	        P1 = v1[ 2 ];
	        P2 = v1[ 3 ];
	        break;
		
		}
	}
	
      
     //T1 = Tnew - Told;
	 //frequency = 1.0 / (Tnew - Told);
     //printf("%f %f %f %f %f %f\n", mu, z, T1, frequency, Vmax, Vmin); 
     //fprintf(fro, "%f %f %f %f %f %f\n", mu, z, T1, frequency, Vmax, Vmin); 
	 

	 if(frequency == 0.0 ) {
		 z_sn = z + 0.0001;
		 z_0 = z_sn + 0.0003;
		 break;
	 }		

	 //printf("%f %f %f %f\n", mu, z, frequency, frequency_old);
	 frequency_old = frequency;
	

	}
	//fclose( fro );

	frequency = frequency_old;
	
    //printf("%f %f %f %f %f %f %f %f\n", d, z_sn, z_0, V, W, P1, P2, frequency);  
    //fprintf(fro, "%f %f %f %f %f %f %f %f\n", d, z_sn, z_0, V, W, P1, P2, frequency);  
    //}

	//fclose( fro );

     //初期条件再設定
	 z = z_sn; //定常カレントの再設定に注意！
     count1 = 0;
      
     t = 0.0;
      
     v[ 0 ] = V;
     v[ 1 ] = W;
     v[ 2 ] = F( 0, t, (double *)v );
     v[ 3 ] = F( 1, t, (double *)v );
     
	 c1[ count1 ][ 0 ] = t;
     c1[ count1 ][ 1 ] = v[ 0 ];
     c1[ count1 ][ 2 ] = v[ 1 ];
     c1[ count1 ][ 3 ] = v[ 2 ];
     c1[ count1 ][ 4 ] = v[ 3 ];
     
	 //printout1( count1 );
	 //printout( t );
	 //printf("%f %f %f %f %f %f %f\n", d, z_sn, v[0], v[1], v[2], v[3], frequency);  
     //printf("%f %f %f %f %f %f %f\n", d, z_sn, c1[0][1], c1[0][2], c1[0][3], c1[0][4], frequency);  
     //printout1( count1 );
	 
	 for( t = H ; t < T1 + H ; t = t + H ) {
		 
		 count1++;
		 c1[ count1 ][ 0 ] = t / T1;
		 runge( t );
		 
		 for(int j = 0 ; j < 4 ; j++) {
			 
			 v[ j ] = v1[ j ];
			 c1[ count1 ][ j + 1 ] = v1[ j ];
		 
		 }

		 n_inf = pow( (1.0 + exp( - (v[ 0 ] - Va) / ka) ), -4.0);
		 b_inf = pow( 1.0 + exp( 10.0 * (v[ 0 ] + 53.3) ) , -1.0/4.0);

		 printf("%f %f %f %f %f\n", t / T1, v[ 0 ], a * v[ 0 ] * v[ 0 ] + b * v[ 0 ] + c,  n_inf + b_inf, 1.0 / (1.0 + exp( - (v[ 0 ] - Va) / ka) ) / (1.0 + exp( - (v[ 0 ] - Va) / ka) )/ (1.0 + exp( - (v[ 0 ] - Va) / ka) )/ (1.0 + exp( - (v[ 0 ] - Va) / ka) ));
		 fprintf(fro, "%f %f %f %f %f\n", t / T1, v[ 0 ], a * v[ 0 ] * v[ 0 ] + b * v[ 0 ] + c, n_inf + b_inf, 1.0 / (1.0 + exp( - (v[ 0 ] - Va) / ka) ) / (1.0 + exp( - (v[ 0 ] - Va) / ka) )/ (1.0 + exp( - (v[ 0 ] - Va) / ka) )/ (1.0 + exp( - (v[ 0 ] - Va) / ka) ));
		 
		 //printout1( count1 );
		 //printout( t );
	 }
	 
	 fclose( fro );

	  //printf("%f %f %f %f %f %f %f %f\n", d, z_sn, z_0, c1[0][1], c1[0][2], c1[0][3], c1[0][4], frequency);  

	 //随伴解導出
      t = 0.0;
      count2 = count1;
	  
	  //初期条件再々設定
      q[ 0 ] = c1[ count2 ][ 1 ];
      q[ 1 ] = c1[ count2 ][ 2 ];
      q[ 2 ] = c1[ count2 ][ 3 ];
      q[ 3 ] = c1[ count2 ][ 4 ];
      
	  //随伴解初期条件設定
	  for(int j = 4 ; j < 6 ; j ++) {
		  //q[ j ]= 0.5; //(double)rand()/((double)(RAND_MAX)+1.0)-0.32;
		  q[ j ] = q[ j - 2 ];
	  } 
 
      //線形解と随伴解の正規直交条件
      CCold = q[ 2 ] * q[ 4 ] + q[ 3 ] * q[ 5 ];
      
      q[ 4 ] = q[ 4 ] / CCold;
      q[ 5 ] = q[ 5 ] / CCold;
      
	  //正規直交化確認
      CC = q[ 2 ] * q[ 4 ] + q[ 3 ] * q[ 5 ];
      //printf("%f\n",CC);
      
      for( t = H ; t < T1 + H ; t = t + H ) {
		  
		  count2--;
		  runge1( t );		  
		  for( int j = 0 ; j < 6 ; j ++ ) {
			  
			  if( 0 <= j && j < 4 ){
				  q[ j ] = c1[ count2 ][ j + 1 ];
			  }	
			  else {
				  q[ j ] = q1[ j ];
				  c1[ count2 ][ j + 1 ] = q[ j ];
			  }
		  }

		  CC = q[ 2 ] * q[ 4 ] + q[ 3 ] * q[ 5 ];
		  c1[ count2 ][ 7 ] = CC;
		  //printf("%f\n", c[ count2 ][ 13 ]);
	  }
	  
	  //全解確認
	  //for( i2 = 0 ; i2 < count1 ; i2 = i2 + 1) {
	  //  //
	 // 	  //printout2( i2 );
	//	  printf("%f %f %f\n", c1[ i2 ][ 0 ], c1[ i2 ][ 1 ], a * c1[ i2 ][ 1 ] * c1[ i2 ][ 1 ] + b * c1[ i2 ][ 1 ] + c ); //, 1.0 / (1.0 + exp ( - (c[ i2 ][ 1 ] - Va) / ka) ) );
	//	  fprintf(fro, "%f %f %f\n", c1[ i2 ][ 0 ], c1[ i2 ][ 1 ], a * c1[ i2 ][ 1 ] * c1[ i2 ][ 1 ] + b * c1[ i2 ][ 1 ] + c ); //, 1.0 / (1.0 + exp( - (c[ i2 ][ 1 ] - Va) / ka ) ) );
//
//		  //printf("%f %f %f %f %f %f %f\n", c1[ i2 ][ 0 ], c1[ i2 ][ 5 ], c1[ i2 ][ 3 ] * c1[ i2 ][ 5 ], c * c1[ i2 ][ 1 ] * c1[ i2 ][ 5 ], - c * c1[ i2 ][ 1 ] * c1[ i2 ][ 1 ] * c1[ i2 ][ 1 ] / 3.0 * c1[ i2 ][ 5 ], - c * c1[ i2 ][ 2 ] * c1[ i2 ][ 5 ], c * c1[ i2 ][ 5 ] * z );
//		  //fprintf(fro, "%f %f %f %f %f %f %f\n", c1[ i2 ][ 0 ], c1[ i2 ][ 5 ], c1[ i2 ][ 3 ] * c1[ i2 ][ 5 ], c * c1[ i2 ][ 1 ] * c1[ i2 ][ 5 ], - c * c1[ i2 ][ 1 ] * c1[ i2 ][ 1 ] * c1[ i2 ][ 1 ] / 3.0 * c1[ i2 ][ 5 ], - c * c1[ i2 ][ 2 ] * c1[ i2 ][ 5 ], c * c1[ i2 ][ 5 ] * z );
//	 }

//	  fclose( fro );

	  //fprintf(fro4, "\n");
	  //}

	  //printf("%d\n", count1);
	  
	  /*
	  //double ZZ_t = 0.0;
	  double ZZ_w1 = 0.0;
	  double ZZ_v1 = 0.0;
	  
	  double ZF_v1 = 0.0;
	  double ZF_v2 = 0.0;

	  double ZF_w1 = 0.0;
	  double WF_w1 = 0.0;
	  double VF_w1 = 0.0;

	  double Z_wW1 = 0.0;
	  double Z_wV1 = 0.0;

	  double AA_1 = 0.0;
	  double AA_2 = 0.0;

	  double dfdd = 0.0;
	  double dfdb = 0.0;

	  double dfdd_o = 0.0;
	  double dfdd_o1 = 0.0;

	  //double ZF_v1 = 0.0;
	  //double ZF_w1 = 0.0;
	  for( i2 = 1 ; i2 < count1 ; i2++ ) {
			//ZZ_t += c1[ i2 ][ 5 ] * ( c1[ i2 + 1 ][ 0 ] - c1[ i2 ][ 0 ] ) ; /// T1; 
			ZZ_v1 += c1[ i2 ][ 5 ] * ( c1[ i2 + 1 ][ 0 ] - c1[ i2 ][ 0 ] ) ; // * frequency;
			ZZ_w1 += c1[ i2 ][ 6 ] * ( c1[ i2 + 1 ][ 0 ] - c1[ i2 ][ 0 ] ) ; // * frequency;
			Z_wW1 += c1[ i2 ][ 6 ] * c1[ i2 ][ 2 ] * ( c1[ i2 + 1 ][ 0 ] - c1[ i2 ][ 0 ] ) ; // * frequency;
			Z_wV1 += c1[ i2 ][ 6 ] * c1[ i2 ][ 1 ] * ( c1[ i2 + 1 ][ 0 ] - c1[ i2 ][ 0 ] ) ; // * frequency;


			//ZZ_v += c1[ i2 ][ 9 ] * ( c1[ i2 + 1 ][ 1 ] - c1[ i2 ][ 1 ] ) ;
			// ZF_v += c1[ i2 ][ 9 ] * c1[ i2 ][ 5 ] * ( c1[ i2 + 1 ][ 0 ] - c1[ i2 ][ 0 ] );
			ZF_v2 += c1[ i2 ][ 5 ] * ( c1[ i2 ][ 3 ] - c2 * z_sn ) * ( c1[ i2 + 1 ][ 0 ] - c1[ i2 ][ 0 ] );
			//ZF_v1 += c1[ i2 ][ 5 ] * c1[ i2 ][ 3 ] * ( c1[ i2 + 1 ][ 0 ] - c1[ i2 ][ 0 ] );
			
			VF_w1 += c1[ i2 ][ 1 ] * c1[ i2 ][ 4 ] * ( c1[ i2 + 1 ][ 0 ] - c1[ i2 ][ 0 ] );
			ZF_v1 += c1[ i2 ][ 5 ] / c1[ i2 ][ 3 ] * ( c1[ i2 + 1 ][ 0 ] - c1[ i2 ][ 0 ] ) ; // * frequency;
			
			WF_w1 += c1[ i2 ][ 2 ] / c1[ i2 ][ 4 ] * ( c1[ i2 + 1 ][ 0 ] - c1[ i2 ][ 0 ] );
			ZF_w1 += c1[ i2 ][ 6 ] * c1[ i2 ][ 4 ] * ( c1[ i2 + 1 ][ 0 ] - c1[ i2 ][ 0 ] ) ; // * frequency;
	  }

	  AA_1 = Z_wW1 - WF_w1 * ZF_w1;
	  AA_2 = Z_wV1 - VF_w1 * ZF_w1;

	  dfdd = - frequency * mu / c2 * ( AA_1 + WF_w1 * ( 1.0 - ZF_v2 - ZZ_v1 * ( c2 * z_sn ) ) );
      dfdb = frequency * mu / c2 * ( AA_2 + VF_w1 * ( 1.0 - ZF_v2 - ZZ_v1 * ( c2 * z_sn ) ) );

	  dfdd_o = - mu / c2 * frequency * Z_wW1;
      dfdd_o1 = (frequency - frequency1)/0.001; 

	  printf("%f %f %f %f %f %f %f %f %f %f %f %f %f\n", d, mu, c2 * z_sn, frequency * 1000.0, ZZ_v1, ZZ_v1 * frequency, ZF_w1, ZF_w1 * frequency, ZF_v2, AA_1, dfdd, dfdd_o, dfdd_o1);
	  fprintf(fro, "%f %f %f %f %f %f %f %f %f %f %f %f %f\n", d, mu, c2 * z_sn, frequency * 1000.0, ZZ_v1, ZZ_v1 * frequency, ZF_w1, ZF_w1 * frequency, ZF_v2, AA_1, dfdd, dfdd_o, dfdd_o1);
	  

	  //printf("%f %f %f %f %f %f %f %f %f\n", mu, z, ZF_v, ZF_v1, frequency * 1000.0, AA_2, AA_3, ZF_w1, ZZ_w);
	  //fprintf( fro, "%f %f %f %f %f %f %f %f\n", mu, z, ZF_v, frequency * 1000.0, AA_2, AA_3, ZF_w1, ZZ_w);
	  
	  //printf("%f %f %f %f %f %f %f %f\n", mu, z, frequency * 1000.0, c * ZZ_v, mu * ZZ_w / c, ZF_v, ZF_w / mu, AA2);
	  //fprintf( fro, "%f %f %f %f %f %f %f %f\n", mu, z, frequency * 1000.0, c * ZZ_v, mu * ZZ_w / c, ZF_v, ZF_w / mu, AA2);
	  
	  //printf("%f %f %f %f %f %f %f %f\n", mu, z, frequency * 1000.0, fabs(frequency - frequency_old) / delta_z, c * ZZ_v, mu * ZZ_w / c, ZF_v, ZF_w / mu);
	  //fprintf( fro, "%f %f %f %f %f %f %f %f\n", mu, z, frequency * 1000.0, fabs(frequency - frequency_old) / delta_z, c * ZZ_v, mu * ZZ_w / c, ZF_v, ZF_w / mu);
	  
	  frequency1 = frequency_old;
	  */

	  //}
	  
	  //fclose( fro );

  //}
	  //位相応答曲線内の面積比導出．
	  //double ZZ_plus = 0.0;
	  //double ZZ_minus = 0.0;
	  //double ZZ = 0.0;
	  //int ip = 0;
	  //int im = 0;
	  //int ip1 = 0;
	  //for( i2 = 1 ; i2 <= count1 - 1 ; i2 = i2 + 1) {
	//	  if( c[ i2 ][ 10 ] * ( 1000.0 / T1 ) > 0 ) {
	//		  ip ++;
	//		  ZZ_plus += c[ i2 ][ 10 ] * ( 1000.0 / T1 ) ;
	//	  }
	//	  else {
	//		  im ++;
	//		  ZZ_minus += c[ i2 ][ 10 ] * ( 1000.0 / T1 ) ;
	//	  }
	//	  ip1 ++;
    //      ZZ += c[ i2 ][ 10 ] * ( 1000.0 / T1 );

	 // }
      

	  //printf( "%f %f %f %f %f %f\n", I_app, mu, ZZ/C/(double)(ip1), frequency, (frequency - frequency1)/0.001, fabs(ZZ_minus) / ZZ_plus );      //( - fabs( ZZ_minus ) + ZZ_plus ) / (double)(im + ip) / C );
	  //fprintf( fro7, "%f %f %f %f %f %f\n", I_app, mu, ZZ/C/(double)(ip1), frequency, (frequency - frequency1)/0.001, fabs(ZZ_minus) / ZZ_plus); //( - fabs( ZZ_minus ) + ZZ_plus ) / (double)(im + ip) / C );
	  //frequency1 = frequency;
	  //fprintf(fro4, "%f %f\n", I_app, ZZ_plus / (double)ip );
	  //fprintf(fro5, "%f %f\n", I_app, fabs( ZZ_minus )  / (double)im );

	  //fprintf(fro4, "%f %f %f\n", I_app, mu, ZZ_plus  );
	  //fprintf(fro5, "%f %f\n", I_app, fabs( ZZ_minus ) );
    
	  //fprintf(fro4, "%f %f %f\n", I_app, ZZ_plus  / ((double)ip + (double)im) );
	  //fprintf(fro5, "%f %f %f\n", I_app, fabs( ZZ_minus ) / ((double)ip + (double)im));

    //}
	  //printf("\n");
	  //fprintf(fro4, "\n");
	  //fprintf(fro5, "\n");
	//}
   //fprintf(fro4,"%f %f %f %f %f %f %f %f %f\n", c[i2][1]/T1, c[i2][2], c[i2][3], c[i2][4], c[i2][5], c[i2][6], c[i2][10], c[i2][14], c[i2][6]*c[i2][10] );
   //fclose(fro4);
   //fclose(fro7);
	  
	  //printf("%d\n", count1);


}

