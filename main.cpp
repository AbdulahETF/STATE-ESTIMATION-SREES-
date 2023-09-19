#include <iostream>
#include <complex>
#include <cmath>

#include <string>
#include <vector>
#include <fstream>
#include <algorithm>

#define N 9
using namespace std;
/*

1. Prođemo naš zadatak i uporedimo sa xml
2. Završiti program

*/


 

 
// Function to get adjoint of A[N][N] in adj[N][N].

 
// Function to calculate and store inverse, returns false if
// matrix is singular

vector<vector<double>> mnozenje_matrica(vector<vector<double>> M1, vector<vector<double>> M2)
{
	int i, j, k;
    
    vector<vector<double>> M;
	
	for(i = 0; i < M1.size(); ++i)
	{
	    	vector<double> V;

		for(j = 0; j < M2[i].size(); ++j)
		{
			V.push_back(0.0);
		}
		
		  M.push_back(V);
	}

	for(i = 0; i < M1.size(); ++i)
	{
		for(j = 0; j < M2[i].size(); ++j)
		{
			for(k=0; k<M1[i].size(); ++k)
			{
				M[i][j] += M1[i][k] * M2[k][j];
			}
		//	if(M[i][j] < 1e-10) {
		//	    M[i][j] = 0.0;
		//	}
		}
	}
	
	
	
	return M;
}


vector<vector<double>> transponuj (vector<vector<double>> M1 ) {
    
    vector<double> V (M1.size(),0);
    vector<vector<double>> T(M1[0].size(),V);

    
	
		for(int i = 0; i < M1.size(); ++i)  {
		    
		for(int j = 0; j < M1[i].size(); ++j)   {
		    
		    T[j][i] = M1[i][j];
		    
		}
		
		}
	return T;
    
}

void kofaktor(const vector<vector<double>> &M,vector<vector<double>> &M_help, int p, int q,int n)
{
    int i = 0, j = 0;
 
    // Looping for each element of the matrix
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            //  Copying into temporary matrix only those
            //  element which are not in given row and
            //  column
            if (row != p && col != q) {
                M_help[i][j++] = M[row][col];
 
                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

double determinant(const vector<vector<double>> &M,int n)
{
   double D = 0; // Initialize result
 
    //  Base case : if matrix contains single element
    if (n == 1)
        return M[0][0];
   vector<double> V (M.size(),0);
    vector<vector<double>> temp(M[0].size(),V);
    
 
    int sign = 1; // To store sign multiplier
 
    // Iterate for each element of first row
    for (int f = 0; f < n; f++) {
        // Getting Cofactor of A[0][f]
        kofaktor(M, temp, 0, f, n);
        D += sign * M[0][f] * determinant(temp, n - 1);
 
        // terms are to be added with alternate sign
        sign = -sign;
    }

    return D;
}


void adjungovana(const vector<vector<double>> &M, vector<vector<double>> &adj)
{
   
  
   
    // temp is used to store cofactors of A[][]
    int sign = 1;
   vector<double> V(M[0].size(),0);
   vector<vector<double>> M_help(M.size(),V);
 
   if (M.size() == 1) {
        adj[0][0] = 1;
        return;
    }
 

 
    for (int i = 0; i < M.size(); i++) {
        for (int j = 0; j < M.size(); j++) {
            // Get cofactor of A[i][j]
            kofaktor(M, M_help, i, j, M.size());
 
            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i + j) % 2 == 0) ? 1 : -1;
 
            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j][i] = (sign) * (determinant(M_help, 8));
        }
    }
}

 

bool inverzna(const vector<vector<double>> &M,vector<vector<double>> &inv) {
    
    // Find determinant of A[][]
    double det = determinant(M, M.size());
    if (det == 0) {
        cout << "Singular matrix, can't find its inverse";
        return false;
    }
 
   vector<double> V(M[0].size(),0.0);
   vector<vector<double>> adj(M.size(),V);
 
    // Find adjoint
    
    adjungovana(M, adj);
 
    // Find Inverse using formula "inverse(A) =
    // adj(A)/det(A)"
    for (int i = 0; i < M.size(); i++)
        for (int j = 0; j < M[0].size(); j++)
            inv[i][j] = adj[i][j] / det;
 
    return true;
    
    
}



int main()
{
    complex<double> Z_s (0.01, 0.04); // serijska impedansa linije 
    complex<double> Y_c (0.0, 0.5); // admitansa kondenzatora u čvoru dva 
    
    //vector<double> z; // vektor mjerenja 
    vector<double> z = {3.58,2.00,1.05,1.92,0.99,-1.62,-0.88,0.12,-0.02,0.997,-2.00,-0.80}; // vektor mjerenja 
   
    double unos; // varijabla za unos 
    vector<string> imena = {"P1","Q1","V1","P13","Q13","P41","Q41","P42","Q42","V2","P3","Q3"};
    
    
      vector<double> rez_vek;
    
   /* for(int i = 0;i<12;i++) {
        
       cout<<imena[i]<<" : ";
       cin>>unos;
       z.push_back(unos);
        
        
    }*/
    
     double delta1 = 0.0;
     double delta2 = 0.0;
     double delta3 = 0.0;
     double delta4 = 0.0;
     
     double V1 = 1.0;
     double V2 = 1.0;
     double V3 = 1.0;
     double V4 = 1.0;
     double epsilon = 1e-3;
     vector<double> X = {delta2,delta3,delta4,V1,V2,V3,V4};
     
     vector<double> unosi = {0.83,0.83,0.417,0.67,0.67,0.67,0.67,1,1,0.417,0.83,0.83};
     int counter = 0;
     // težinski vektor 
     cout<<"Unesite vrijednosti težinskog vektora: "<<endl;
     vector<vector<double>> W;
     vector<double> v11; // prazan vektor
     for(int i = 0;i<12;i++) {
         W.push_back(v11);
        for(int j = 0;j<12;j++) {
         
         if(i == j) {
         cout<<i+1<<". vrijednost : ";
        // cin>>unos;
         //W[i].push_back(unos);
         W[i].push_back(unosi[counter]);
         counter++;
         }
         else {
             W[i].push_back(0.0);
         }
        }
     }
     
     cout<<endl<<"W: "<<endl;
     for(int i = 0;i<12;i++) {
        for(int j = 0;j<12;j++) {
         
         cout<<W[i][j]<<" ";
         
        }
        cout<<endl;
     }
     
    vector<vector<complex<double>>> Y; // matrica admitansi sistema
    
    
    for(int i = 0;i<4;i++) {
    vector<complex<double>> VC; 
    
    for(int j = 0;j<4;j++) {
        
        if ( i == j && i != 1) {
            
                       VC.push_back(2.0/Z_s);

        }
        else if( i==j && i==1) {
            VC.push_back((2.0/Z_s) + Y_c);

        }
        else if(i == 0 && j == 2 || i == 0 && j == 3 || i == 1 && j == 2 || i == 1 && j==3 || i == 2 && j == 0 || i == 2 && j == 1 || i == 3 && j==0 || i == 3 && j==1){
            
            VC.push_back(-1.0/Z_s);
        }
        else {
            VC.push_back(complex<double>(0,0));

            
        }
        
        
    
        
    }
       Y.push_back(VC);
    }


    for(int i = 0;i<4;i++) {
        
        for(int j = 0;j<4;j++) {
            
            cout<<Y[i][j]<<" ";
        }
        cout<<endl;
        
    }
    int counter_iteracija = 0;
    for(int z_j=0;z_j<5;z_j++) {
    vector<vector<double>> C;
    for(int i = 0; i<2;i++) {
    vector<double> pomocni;
    if( i == 0) {
    double x1 = -1.0*(abs(Y[1][1])*abs(V2)*sin(arg(Y[1][1])+delta2));
    
    pomocni.push_back(x1);
    double x2 = -1*(abs(Y[1][2])*abs(V3)*sin(arg(Y[1][2])+delta3));
    
    pomocni.push_back(x2);
    
    double x3 = -1*(abs(Y[1][3])*abs(V4)*sin(arg(Y[1][3])+delta4));
   
    pomocni.push_back(x3);
    pomocni.push_back(0.0);
    
    
    double x4 = -1*(abs(Y[1][1])*cos(arg(Y[1][1]+delta2)));
   
   
    pomocni.push_back(x4);    
    
    double x5 = -1*(abs(Y[1][2])*cos(arg(Y[1][2]+delta3)));
     
    
    
    pomocni.push_back(x5);    


    double x6 = -1*(abs(Y[1][3])*cos(arg(Y[1][3]+delta4)));
    
     
    pomocni.push_back(x6);
    }
    
    else {
        
    double x1 = -1*(abs(Y[1][1])*abs(V2)*cos(arg(Y[1][1])+delta2));
     
    pomocni.push_back(x1);
    double x2 = -1*(abs(Y[1][2])*abs(V3)*cos(arg(Y[1][2])+delta3));
  
    pomocni.push_back(x2);
    
    double x3 = -1*(abs(Y[1][3])*abs(V4)*cos(arg(Y[1][3])+delta4));
    
    
    pomocni.push_back(x3);
    pomocni.push_back(0.0);
    
    
    double x4 = -1*(abs(Y[1][1])*sin(arg(Y[1][1]+delta2)));
    
     
    pomocni.push_back(x4);    
    
    double x5 = -1*(abs(Y[1][2])*sin(arg(Y[1][2]+delta3)));
     
    pomocni.push_back(x5);    


    double x6 = -1*(abs(Y[1][3])*sin(arg(Y[1][3]+delta4)));
    
    pomocni.push_back(x6);
        
        
    }
    

    
     C.push_back(pomocni);
    
    
    
    }
    
    
    cout<<endl<<"C: "<<endl;
     for(int i = 0;i<2;i++) {
        for(int j = 0;j<7;j++) {
         
         cout<<C[i][j]<<" ";
         
        }
        cout<<endl;
     }
     
     
     complex<double> y(5.88235,-23.52941);
      complex<double> y_neg(-5.88235,23.52941);
    
     // y mora biti matrica??? kakva nemamo pojma
    
    
    
    double P1 = abs(V1)*abs(V1)*abs(Y[0][0])*cos(arg(Y[0][0])) + abs(V1)*abs(V2)*abs(Y[0][1])*cos(arg(Y[0][1])-delta1+delta2)+abs(V1)*abs(V3)*abs(Y[0][2])*cos(arg(Y[0][2]-delta1+delta3)) + abs(V1)*abs(V4)*abs(Y[0][3])*cos(arg(Y[0][3]-delta1+delta4)) ;
    
    double Q1 =  -1*pow(abs(V1),2)*abs(Y[0][0])*sin(arg(Y[0][0])) - abs(V1)*abs(V2)*abs(Y[0][1])*sin(arg(Y[0][1])-delta1+delta2)-abs(V1)*abs(V3)*abs(Y[0][2])*sin(arg(Y[0][2]-delta1+delta3)) - abs(V1)*abs(V4)*abs(Y[0][3])*sin(arg(Y[0][3]-delta1+delta4)); 
    
    double P3 = pow(abs(V3),2)*abs(Y[2][2])*cos(arg(Y[2][2])) + abs(V3)*abs(V1)*abs(Y[2][0])*cos(arg(Y[2][0])-delta3+delta1)+abs(V3)*abs(V2)*abs(Y[2][1])*cos(arg(Y[2][1]-delta3+delta2)) + abs(V3)*abs(V4)*abs(Y[2][3])*cos(arg(Y[2][3]-delta3+delta4)); 
    
    double Q3 =  -1*pow(abs(V3),2)*abs(Y[2][2])*sin(arg(Y[2][2])) - abs(V3)*abs(V1)*abs(Y[2][0])*sin(arg(Y[2][0])-delta3+delta1)-abs(V3)*abs(V2)*abs(Y[2][1])*sin(arg(Y[2][1]-delta3+delta2)) - abs(V3)*abs(V4)*abs(Y[2][3])*sin(arg(Y[2][3]-delta3+delta4)); 
    
    
    double P13 = pow(abs(V1),2)*abs(y)*cos(arg(y)) + abs(V1)*abs(V3)*abs(y_neg)*cos(arg(y_neg)-delta1+delta3);
    
    //double P13 = pow(abs(V1),2)*abs(Y[0][2])*cos(arg(Y[0][2])) + abs(V1)*abs(V3)*abs(arg(complex<double>(-1,0)*Y[0][2]))*cos(arg(Y[0][2])-delta1+delta3);
    
    
    double Q13 = -1*pow(abs(V1),2)*abs(y)*sin(arg(y)) - abs(V1)*abs(V3)*abs(y_neg)*sin(arg(y_neg)-delta1+delta3);
  //  double Q13 = -1*pow(abs(V1),2)*abs(Y[0][2])*sin(arg(Y[0][2])) - abs(V1)*abs(V3)*abs(complex<double>(-1,0)*Y[0][2])*sin(arg(Y[0][2])-delta1+delta3);

    
    double P41 = pow(abs(V4),2)*abs(y)*cos(arg(y)) + abs(V1)*abs(V4)*abs(y_neg)*cos(arg(y_neg)-delta4+delta1);
   // double P41 = pow(abs(V4),2)*abs(Y[3][0])*cos(arg(Y[3][0])) + abs(V1)*abs(V4)*abs(complex<double>(-1,0)*Y[3][0])*cos(arg(Y[3][0])-delta4+delta1);




    double Q41 = -1*pow(abs(V4),2)*abs(y)*sin(arg(y)) - abs(V1)*abs(V4)*abs(y_neg)*sin(arg(y_neg)-delta4+delta1);

   // double Q41 = -1*pow(abs(V4),2)*abs(Y[3][0])*sin(arg(Y[3][0])) - abs(V1)*abs(V4)*abs(complex<double>(-1,0)*Y[3][0])*sin(arg(Y[3][0])-delta4+delta1);
   

        double P42 = pow(abs(V4),2)*abs(y)*cos(arg(y)) + abs(V2)*abs(V4)*abs(y_neg)*cos(arg(y_neg)-delta4+delta2);

    
   // double P42 = pow(abs(V4),2)*abs(Y[3][1])*cos(arg(Y[3][1])) + abs(V2)*abs(V4)*abs(complex<double>(-1,0)*Y[3][1])*cos(arg(Y[3][1])-delta4+delta2);
        
        



       double Q42 = -1*pow(abs(V4),2)*abs(y)*sin(arg(y)) - abs(V2)*abs(V4)*abs(y_neg)*sin(arg(y_neg)-delta4+delta2);

   // double Q42 = -1*pow(abs(V4),2)*abs(Y[3][1])*sin(arg(Y[3][1])) - abs(V2)*abs(V4)*abs(complex<double>(-1,0)*Y[3][1])*sin(arg(Y[3][1])-delta4+delta2);
    

 
    
    
    vector<double> h;
    
    h.push_back(P1);
    h.push_back(Q1);
    h.push_back(V1);
    h.push_back(P13);

    h.push_back(Q13);
    h.push_back(P41);
    h.push_back(Q41);
    h.push_back(P42);

    h.push_back(Q42);
    h.push_back(V2);
    h.push_back(P3);
    h.push_back(Q3);
    
     cout<<endl<<endl;
    cout<<"h:"<<endl<<endl;
    for(int i = 0;i<h.size();i++) {
        
       
            cout<<h[i]<<" ";
            
        
        
        
        
    }
    cout<<endl;
    vector<vector<double>> H;
    
    vector<double> prvi_red_H;
    vector<double> drugi_red_H;
    vector<double> treci_red_H;
    vector<double> cetvrti_red_H;
    vector<double> peti_red_H;
    vector<double> sesti_red_H;
    vector<double> sedmi_red_H;
    vector<double> osmi_red_H;
    vector<double> deveti_red_H;
    vector<double> deseti_red_H;
    vector<double> jedanaesti_red_H;
    vector<double> dvanaesti_red_H;
    
    
    
    double dP1delta2 = -1.0*abs(V1)*abs(V2)*abs(Y[0][1])*sin(arg(Y[0][1]) - delta1 + delta2);
    double dP1delta3 = -1*abs(V1)*abs(V3)*abs(Y[0][2])*sin(arg(Y[0][2]) - delta1 + delta3);
    double dP1delta4 = -1*abs(V1)*abs(V4)*abs(Y[0][3])*sin(arg(Y[0][3]) - delta1 + delta4);
  
    
    prvi_red_H.push_back(dP1delta2);
    prvi_red_H.push_back(dP1delta3);
    prvi_red_H.push_back(dP1delta4);
    
    double dP1V1 = 2*abs(V1)*abs(Y[0][0])*cos(arg(Y[0][0])) + abs(V2)*abs(Y[0][1])*cos(arg(Y[0][1])-delta1 + delta2) + abs(V3)*abs(Y[0][2])*cos(arg(Y[0][2])-delta1 + delta3) + abs(V4)*abs(Y[0][3])*cos(arg(Y[0][3] - delta1 + delta4));
    double dP1V2 = abs(V1)*abs(Y[0][1])*cos(arg(Y[0][1]) - delta1 + delta2);
    double dP1V3 = abs(V1)*abs(Y[0][2])*cos(arg(Y[0][2]) - delta1 + delta3);
    double dP1V4 = abs(V1)*abs(Y[0][3])*cos(arg(Y[0][3]) - delta1 + delta4);
    
    
    prvi_red_H.push_back(dP1V1);
    prvi_red_H.push_back(dP1V2);
    prvi_red_H.push_back(dP1V3);
    prvi_red_H.push_back(dP1V4);
    
    
    double dQ1delta2 = -1*abs(V1)*abs(V2)*abs(Y[0][1])*cos(arg(Y[0][1]) - delta1 + delta2);
    double dQ1delta3 = -1*abs(V1)*abs(V3)*abs(Y[0][2])*cos(arg(Y[0][2]) - delta1 + delta3);
    double dQ1delta4 = -1*abs(V1)*abs(V4)*abs(Y[0][3])*cos(arg(Y[0][3]) - delta1 + delta4);
    
    drugi_red_H.push_back(dQ1delta2);
    drugi_red_H.push_back(dQ1delta3);
    drugi_red_H.push_back(dQ1delta4);

    
    double dQ1V1 =  -2*abs(V1)*abs(Y[0][0])*sin(arg(Y[0][0])) - abs(V2)*abs(Y[0][1])*sin(arg(Y[0][1])-delta1 + delta2) - abs(V3)*abs(Y[0][2])*sin(arg(Y[0][2])-delta1 + delta3) - abs(V4)*abs(Y[0][3])*sin(arg(Y[0][3] - delta1 + delta4));
    double dQ1V2 =   -1*abs(V1)*abs(Y[0][1])*sin(arg(Y[0][1]) - delta1 + delta2); 
    double dQ1V3 =   -1*abs(V1)*abs(Y[0][2])*sin(arg(Y[0][2]) - delta1 + delta3);
    double dQ1V4 =   -1*abs(V1)*abs(Y[0][3])*sin(arg(Y[0][3]) - delta1 + delta4);
    
    
    drugi_red_H.push_back(dQ1V1);
    drugi_red_H.push_back(dQ1V2);
    drugi_red_H.push_back(dQ1V3);
    drugi_red_H.push_back(dQ1V4);

    
    treci_red_H.push_back(0);
    treci_red_H.push_back(0);
    treci_red_H.push_back(0);
    treci_red_H.push_back(1);
    treci_red_H.push_back(0);
    treci_red_H.push_back(0);
    treci_red_H.push_back(0);

    
    double dP13delta2 = 0;
    double dP13delta3 = -1*abs(V1)*abs(V3)*abs(y_neg)*sin(arg(y_neg) - delta1 + delta3);
   
     //  double dP13delta3 = -1*abs(V1)*abs(V3)*abs(complex<double>(-1,0)*Y[0][2])*sin(arg(complex<double>(-1,0)*Y[0][2]) - delta1 + delta3);

    double dP13delta4 = 0;
    
/*    if(abs(dP13delta3) < 1e-10) {
        dP13delta3 = 0;
    }
  */  
    cetvrti_red_H.push_back(dP13delta2);
    cetvrti_red_H.push_back(dP13delta3);
    cetvrti_red_H.push_back(dP13delta4);
    
    double dP13V1 = 2*abs(V1)*abs(y)*cos(arg(y)) + abs(V3)*abs(y_neg)*cos(arg(y_neg) - delta1 + delta3);

   // double dP13V1 = 2*abs(V1)*abs(Y[0][2])*cos(arg(Y[0][2])) + abs(V3)*abs(complex<double>(-1,0)*Y[0][2])*cos(arg(complex<double>(-1,0)*Y[0][2]) - delta1 + delta3);
    double dP13V2 = 0;
  double dP13V3 = abs(V1)*abs(y_neg)*cos(arg(y_neg)-delta1 + delta3);

    //double dP13V3 = abs(V1)*abs(complex<double>(-1,0)*Y[0][2])*cos(arg(complex<double>(-1,0)*Y[0][2])-delta1 + delta3);
    double dP13V4 = 0;
    
   /* if(abs(dP13V1) < 1e-10) {
        dP13V1 = 0;
    }
    
    if(abs(dP13V3) < 1e-10) {
        dP13V3 = 0;
    }*/
    
    cetvrti_red_H.push_back(dP13V1);
    cetvrti_red_H.push_back(dP13V2);
    cetvrti_red_H.push_back(dP13V3);
    cetvrti_red_H.push_back(dP13V4);
    
    
    double dQ13delta2 = 0;
    double dQ13delta3 = -1*abs(V1)*abs(V3)*abs(y_neg)*cos(arg(y_neg) - delta1 + delta3);
    
   // double dQ13delta3 = -1*abs(V1)*abs(V3)*abs(complex<double>(-1,0)*Y[0][2])*cos(arg(complex<double>(-1,0)*Y[0][2]) - delta1 + delta3);
    double dQ13delta4 = 0;
    
  /*  if(abs(dQ13delta3) < 1e-10) {
        dQ13delta3 = 0;
    }*/
    
    peti_red_H.push_back(dQ13delta2);
    peti_red_H.push_back(dQ13delta3);
    peti_red_H.push_back(dQ13delta4);
       double dQ13V1 = -2*abs(V1)*abs(y)*sin(arg(y)) - abs(V3)*abs(y_neg)*sin(arg(y_neg) - delta1 + delta3); 

    //double dQ13V1 = -2*abs(V1)*abs(Y[0][2])*sin(arg(Y[0][2])) - abs(V3)*abs(complex<double>(-1,0)*Y[0][2])*sin(arg(complex<double>(-1,0)*Y[0][2]) - delta1 + delta3); 
    double dQ13V2 = 0;
    
        double dQ13V3 = -1*abs(V1)*abs(y_neg)*sin(arg(y_neg)-delta1 + delta3);

    //double dQ13V3 = -1*abs(V1)*abs(complex<double>(-1,0)*Y[0][2])*sin(arg(complex<double>(-1,0)*Y[0][2])-delta1 + delta3);
    double dQ13V4 = 0;
    
   /* if(abs(dQ13V1) < 1e-10) {
        dQ13V1 = 0;
    }
    if(abs(dQ13V3) < 1e-10) {
        dQ13V3 = 0;
    }*/
    
    peti_red_H.push_back(dQ13V1);
    peti_red_H.push_back(dQ13V2);
    peti_red_H.push_back(dQ13V3);
    peti_red_H.push_back(dQ13V4);
    
    
    double dP41delta2 = 0;
    
    
    double dP41delta3 = 0;
        double dP41delta4 = abs(V1)*abs(V4)*abs(y_neg)*sin(arg(y_neg)-delta4 + delta1);

   // double dP41delta4 = abs(V1)*abs(V4)*abs(complex<double>(-1,0)*Y[3][0])*sin(arg(complex<double>(-1,0)*Y[3][0])-delta4 + delta1);
    
   /* if(abs(dP41delta4) < 1e-10) {
        dP41delta4 = 0;
    }*/
    
    sesti_red_H.push_back(dP41delta2);
    sesti_red_H.push_back(dP41delta3);
    sesti_red_H.push_back(dP41delta4);
        double dP41V1 = abs(V4)*abs(y_neg)*cos(arg(y_neg)-delta4 + delta1);

//    double dP41V1 = abs(V4)*abs(complex<double>(-1,0)*Y[3][0])*cos(arg(complex<double>(-1,0)*Y[3][0])-delta4 + delta1);
    double dP41V2 = 0;
    double dP41V3 = 0;
    
    double dP41V4 = 2*abs(V4)*abs(y)*cos(arg(y)) + abs(V1)*abs(y_neg)*cos(arg(y_neg)-delta4 + delta1);
   
   //double dP41V4 = 2*abs(V4)*abs(Y[3][0])*cos(arg(Y[3][0])) + abs(V1)*abs(complex<double>(-1,0)*Y[3][0])*cos(arg(complex<double>(-1,0)*Y[3][0])-delta4 + delta1);
    
  /*  if(abs(dP41V1) < 1e-10) {
        dP41V1 = 0;
    }
    if(abs(dP41V4) <1e-10) {
        dP41V4 = 0;
    }*/
    
    sesti_red_H.push_back(dP41V1);
    sesti_red_H.push_back(dP41V2);
    sesti_red_H.push_back(dP41V3);
    sesti_red_H.push_back(dP41V4);
    
    double dQ41delta2 = 0;
    double dQ41delta3 = 0; 
    double dQ41delta4 = abs(V1)*abs(V4)*abs(y_neg)*cos(arg(y_neg)-delta4 + delta1);
    //double dQ41delta4 = abs(V1)*abs(V4)*abs(complex<double>(-1,0)*Y[3][0])*cos(arg(complex<double>(-1,0)*Y[3][0])-delta4 + delta1);
  
   /* if(abs(dQ41delta4) < 1e-10) {
        dQ41delta4 = 0;
    }*/
  
    sedmi_red_H.push_back(dQ41delta2);
    sedmi_red_H.push_back(dQ41delta3);
    sedmi_red_H.push_back(dQ41delta4);
    
    double dQ41V1 =  -1*abs(V4)*abs(y_neg)*sin(arg(y_neg)-delta4 + delta1);
    
    //double dQ41V1 =  -1*abs(V4)*abs(complex<double>(-1,0)*Y[3][0])*sin(arg(complex<double>(-1,0)*Y[3][0])-delta4 + delta1);
    double dQ41V2 = 0;
    double dQ41V3 = 0;
    double dQ41V4 = -2*abs(V4)*abs(y)*sin(arg(y)) - abs(V1)*abs(y_neg)*sin(arg(y_neg)-delta4 + delta1);
    //double dQ41V4 = -2*abs(V4)*abs(y)*sin(arg(Y[3][0])) - abs(V1)*abs(complex<double>(-1,0)*Y[3][0])*sin(arg(complex<double>(-1,0)*Y[3][0])-delta4 + delta1);
    
  /*  if(abs(dQ41V1) < 1e-10) {
        dQ41V1 = 0;
    }
    if(abs(dQ41V4) < 1e-10) {
        dQ41V4 = 0;
    }*/
    sedmi_red_H.push_back(dQ41V1);
    sedmi_red_H.push_back(dQ41V2);
    sedmi_red_H.push_back(dQ41V3);
    sedmi_red_H.push_back(dQ41V4);
        
 
     
    double dP42delta2 =  -1*abs(V2)*abs(V4)*abs(y_neg)*sin(arg(y_neg)-delta4 + delta2);
    
//    double dP42delta2 =  -1*abs(V2)*abs(V4)*abs(complex<double>(-1,0)*Y[3][1])*sin(arg(complex<double>(-1,0)*Y[3][1])-delta4 + delta2);
    double dP42delta3 = 0;
    
    double dP42delta4 = abs(V2)*abs(V4)*abs(y_neg)*sin(arg(y_neg)-delta4 + delta2);
    
    //double dP42delta4 = abs(V2)*abs(V4)*abs(complex<double>(-1,0)*Y[3][1])*sin(arg(complex<double>(-1,0)*Y[3][1])-delta4 + delta2);
    
   /* if(abs(dP42delta2) < 1e-10) {
        dP42delta2 = 0;
    }
    if(abs(dP42delta4) < 1e-10) {
        dP42delta4 = 0;
    }*/
    
    osmi_red_H.push_back(dP42delta2);
    osmi_red_H.push_back(dP42delta3);
    osmi_red_H.push_back(dP42delta4);
    
    double dP42V1 = 0;
    double dP42V2 = abs(V4)*abs(y_neg)*cos(arg(y_neg)-delta4 + delta2);
   
  // double dP42V2 = abs(V4)*abs(complex<double>(-1,0)*Y[3][1])*cos(arg(complex<double>(-1,0)*Y[3][1])-delta4 + delta2);
    double dP42V3 = 0; 
    double dP42V4 = 2*abs(V4)*abs(y)*cos(arg(y)) + abs(V2)*abs(y_neg)*cos(arg(y_neg)-delta4 + delta2);
   
   //double dP42V4 = 2*abs(V4)*abs(y)*cos(arg(Y[3][1])) + abs(V2)*abs(complex<double>(-1,0)*Y[3][1])*cos(arg(complex<double>(-1,0)*Y[3][1])-delta4 + delta2);
   /* if(abs(dP42V2) <1e-10) {
        dP42V2 = 0;
    }
    if(abs(dP42V4) < 1e-10) {
        dP42V4 = 0;
    }*/
    osmi_red_H.push_back(dP42V1);
    osmi_red_H.push_back(dP42V2);
    osmi_red_H.push_back(dP42V3);
    osmi_red_H.push_back(dP42V4);
    
    
    double dQ42delta2 = -1*abs(V2)*abs(V4)*abs(y_neg)*cos(arg(y_neg)-delta4 + delta2);
//    double dQ42delta2 = -1*abs(V2)*abs(V4)*abs(complex<double>(-1,0)*Y[3][1])*cos(arg(complex<double>(-1,0)*Y[3][1])-delta4 + delta2);
    double dQ42delta3 = 0;
    
    double dQ42delta4 =  abs(V2)*abs(V4)*abs(y_neg)*cos(arg(y_neg)-delta4 + delta2);
 //   double dQ42delta4 =  abs(V2)*abs(V4)*abs(complex<double>(-1,0)*Y[3][1])*cos(arg(complex<double>(-1,0)*Y[3][1])-delta4 + delta2);
  /*  if(abs(dQ42delta2) < 1e-10) {
        dQ42delta2 = 0;
    }
    if(abs(dQ42delta4) < 1e-10) {
        dQ42delta4 = 0;
    }*/
    deveti_red_H.push_back(dQ42delta2);
    deveti_red_H.push_back(dQ42delta3);
    deveti_red_H.push_back(dQ42delta4);

    
    double dQ42V1 = 0;
    double dQ42V2 = -1*abs(V4)*abs(y_neg)*sin(arg(y_neg)-delta4 + delta2);
  //  double dQ42V2 = -1*abs(V4)*abs(complex<double>(-1,0)*Y[3][1])*sin(arg(complex<double>(-1,0)*Y[3][1])-delta4 + delta2);
    double dQ42V3 = 0;
    double dQ42V4 = -2*abs(V4)*abs(y)*sin(arg(y)) - abs(V2)*abs(y_neg)*sin(arg(y_neg)-delta4 + delta2);
   // double dQ42V4 = -2*abs(V4)*abs(Y[3][1])*sin(arg(Y[3][1])) - abs(V2)*abs(complex<double>(-1,0)*Y[3][1])*sin(arg(complex<double>(-1,0)*Y[3][1])-delta4 + delta2);
  /*  if(abs(dQ42V2) < 1e-10) {
        dQ42V2 = 0;
    }
    if(abs(dQ42V4) <1e-10) {
        dQ42V4 = 0;
    }*/
    deveti_red_H.push_back(dQ42V1);
    deveti_red_H.push_back(dQ42V2);
    deveti_red_H.push_back(dQ42V3);
    deveti_red_H.push_back(dQ42V4);
    
    
    
    deseti_red_H.push_back(0);
    deseti_red_H.push_back(0);
    deseti_red_H.push_back(0);
    deseti_red_H.push_back(0);
    deseti_red_H.push_back(1);
    deseti_red_H.push_back(0);
    deseti_red_H.push_back(0);
    
    
    double dP3delta2 = -1*abs(V3)*abs(V2)*abs(Y[2][1])*sin(arg(Y[2][1]) - delta3 + delta2);
    double dP3delta3 = abs(V1)*abs(V3)*abs(Y[2][0])*sin(arg(Y[2][0]) - delta3 + delta1) + abs(V3)*abs(V2)*abs(Y[2][1])*sin(arg(Y[2][1]) - delta3 + delta2) + abs(V3)*abs(V4)*abs(Y[2][3])*sin(arg(Y[2][3]) - delta3 + delta4);
    double dP3delta4 =  -1*abs(V3)*abs(V4)*abs(Y[2][3])*sin(arg(Y[2][3]) - delta3 + delta4);
    /*if(abs(dP3delta2) < 1e-10) {
        dP3delta2 = 0;
    }
    if(abs(dP3delta3) < 1e-10) {
        dP3delta3 = 0;
    }
    if(abs(dP3delta4) < 1e-10) {
        dP3delta4 = 0;
    }*/
    jedanaesti_red_H.push_back(dP3delta2);
    jedanaesti_red_H.push_back(dP3delta3);
    jedanaesti_red_H.push_back(dP3delta4);
    
    
    double dP3V1 =  abs(V3)*abs(Y[2][0])*cos(arg(Y[2][0]) - delta3 + delta1);
    double dP3V2 =  abs(V3)*abs(Y[2][1])*cos(arg(Y[2][1]) - delta3 + delta2);
    double dP3V3 =  2*abs(V3)*abs(Y[2][2])*cos(arg(Y[2][2])) + abs(V1)*abs(Y[2][0])*cos(arg(Y[2][0]) - delta3 + delta1) + abs(V2)*abs(Y[2][1])*cos(arg(Y[2][1]) - delta3 + delta2) + abs(V4)*abs(Y[2][3])*cos(arg(Y[2][3]) - delta3 + delta4);
    double dP3V4 = abs(V3)*abs(Y[2][3])*cos(arg(Y[2][3]) - delta3 + delta4);
    /*if(abs(dP3V1) < 1e-10) {
        dP3V1 = 0;
    }
    
    if(abs(dP3V2) < 1e-10) {
        dP3V2 = 0;
    }
    
    if(abs(dP3V3) < 1e-10) {
        dP3V3 = 0;
    }
    
    if(abs(dP3V4) < 1e-10) {
        dP3V4 = 0;
    }
    */
    jedanaesti_red_H.push_back(dP3V1);
    jedanaesti_red_H.push_back(dP3V2);
    jedanaesti_red_H.push_back(dP3V3);
    jedanaesti_red_H.push_back(dP3V4);
    
    double dQ3delta2 = -1*abs(V3)*abs(V2)*abs(Y[2][1])*cos(arg(Y[2][1]) - delta3 + delta2);
    double dQ3delta3 = abs(V3)*abs(V1)*abs(Y[2][0])*cos(arg(Y[2][0]) - delta3 + delta1) + abs(V3)*abs(V2)*abs(Y[2][1])*cos(arg(Y[2][1]) - delta3 + delta2) + abs(V3)*abs(V4)*abs(Y[2][3])*cos(arg(Y[2][3]) - delta3 + delta4);
    double dQ3delta4 = -1*abs(V3)*abs(V4)*abs(Y[2][3])*cos(arg(Y[2][3]) - delta3 + delta4);
  /*  if(abs(dQ3delta2) <1e-10) {
        dQ3delta2 = 0;
    }
    
    if(abs(dQ3delta3) <1e-10) {
        dQ3delta3 = 0;
    }
    
    if(abs(dQ3delta4) <1e-10) {
        dQ3delta4 = 0;
    }*/
    dvanaesti_red_H.push_back(dQ3delta2);
    dvanaesti_red_H.push_back(dQ3delta3);
    dvanaesti_red_H.push_back(dQ3delta4);
    
    
    double dQ3V1 = -1*abs(V3)*abs(Y[2][0])*sin(arg(Y[2][0]) - delta3 + delta1);
    double dQ3V2 = -1*abs(V3)*abs(Y[2][1])*sin(arg(Y[2][1]) - delta3 + delta2);
    double dQ3V3 = -2*abs(V3)*abs(Y[2][2])*sin(arg(Y[2][2])) - abs(V1)*abs(Y[2][0])*sin(arg(Y[2][0]) - delta3 + delta1) -  abs(V2)*abs(Y[2][1])*sin(arg(Y[2][1]) - delta3 + delta2) - abs(V4)*abs(Y[2][3])*sin(arg(Y[2][3]) - delta3 + delta4);
    double dQ3V4 = -1*abs(V3)*abs(Y[2][3])*sin(arg(Y[2][3]) - delta3 + delta4);
   /* 
    if(abs(dQ3V1) < 1e-10) {
        dQ3V1 = 0;
    }
    if(abs(dQ3V2) < 1e-10) {
        dQ3V2 = 0;
    }
    if(abs(dQ3V3) < 1e-10) {
        dQ3V3 = 0;
    }
    if(abs(dQ3V4) < 1e-10) {
        dQ3V4 = 0;
    }
    */
    dvanaesti_red_H.push_back(dQ3V1);
    dvanaesti_red_H.push_back(dQ3V2);
    dvanaesti_red_H.push_back(dQ3V3);
    dvanaesti_red_H.push_back(dQ3V4);
    
    
    H.push_back(prvi_red_H);
    H.push_back(drugi_red_H);
    H.push_back(treci_red_H);  
    H.push_back(cetvrti_red_H);
    H.push_back(peti_red_H);
    H.push_back(sesti_red_H);
    H.push_back(sedmi_red_H);
    H.push_back(osmi_red_H);
    H.push_back(deveti_red_H);
    H.push_back(deseti_red_H);
    H.push_back(jedanaesti_red_H);
    H.push_back(dvanaesti_red_H);
    
   /* cout<<endl<<endl;
    cout<<"H:"<<endl<<endl;
    for(int i = 0;i<12;i++) {
        
        for(int j = 0;j<7;j++) {
            
            cout<<H[i][j]<<" ";
            
        }
        
        cout<<endl;
        
    }*/
    
   /* cout<<endl<<endl;
    vector<vector<double>> M1 = { {3,1,2} ,{0,0,1} };
    vector<vector<double>> M2 = { {2,0,0} , {1,0,1}, {3,0,1}};
    
     for(int i=0;i<M1.size();i++) {
        
        for(int j = 0;j<M1[0].size();j++) {
            
            cout<<M1[i][j]<<" ";
            
        }
        cout<<endl;
    }
    
    cout<<endl<<endl;
   */ 
   // vector<vector<double>> M = mnozenje_matrica(M1,M2);
 //   vector<vector<double>> M = transponuj(M1);
    
   /* for(int i=0;i<M.size();i++) {
        
        for(int j = 0;j<M[0].size();j++) {
            
            cout<<M[i][j]<<" ";
            
        }
        cout<<endl;
    }*/
    
    
    // transponovanje H
   
       vector<vector<double>> H_T = transponuj(H);
  cout<<endl<<endl;
    cout<<"H_T:"<<endl<<endl;
    for(int i = 0;i<H_T.size();i++) {
        
        for(int j = 0;j<H_T[i].size();j++) {
            
            cout<<H_T[i][j]<<" ";
            
        }
        
        cout<<endl;
        
    }
    
    
    // množenje H_T sa W sa H
    
    vector<vector<double>> H_T_W = mnozenje_matrica(H_T,W);
   
     cout<<endl<<endl;
    cout<<"H_T_W:"<<endl<<endl;
    for(int i = 0;i<H_T_W.size();i++) {
        
        for(int j = 0;j<H_T_W[i].size();j++) {
            
            cout<<H_T_W[i][j]<<" ";
            
        }
        
        cout<<endl;
        
    }
    
  vector<vector<double>> H_T_W_H = mnozenje_matrica(H_T_W,H);
  
  
  cout<<endl<<endl;
    cout<<"H_T_W_H:"<<endl<<endl;
    for(int i = 0;i<H_T_W_H.size();i++) {
        
        for(int j = 0;j<H_T_W_H[i].size();j++) {
           
            cout<<H_T_W_H[i][j]<<" ";
            
        }
        
        cout<<endl;
        
    }
    
    
  
    
   
  vector<vector<double>> C_T;
  C_T = transponuj(C);
  
  
  vector<vector<double>>M_M;
  
  
  int brojac_C_T_redovi = 0;
  int brojac_C_T_kolone = 0;
  
  int brojac_C_redovi = 0;
  int brojac_C_kolone = 0;
  
 

  for(int i=0;i<H_T_W_H.size()+C.size();i++) {
      
      vector<double> V_help;
      
      for(int j=0;j<H_T_W_H[0].size()+C_T[0].size();j++) {
          
          
          if(i<H_T_W_H.size() && j<H_T_W_H[0].size()) {
              
              V_help.push_back(H_T_W_H[i][j]);
              continue;
          }
          
          else if(i<H_T_W_H.size() && j>=H_T_W_H[0].size()) {
              
              V_help.push_back(C_T[brojac_C_T_redovi][brojac_C_T_kolone]);
              brojac_C_T_kolone++;
              continue;
          }
          
          else if(i>=H_T_W_H.size() && j<H_T_W_H[0].size()) {
                        
              V_help.push_back(C[brojac_C_redovi][brojac_C_kolone]);
              brojac_C_kolone++;
              continue;
          }
          else {
              
              V_help.push_back(0.0);
              continue;
          }
          
          
      }
                M_M.push_back(V_help);

      brojac_C_T_redovi++;
      brojac_C_T_kolone = 0;
      if(i>=H_T_W_H.size()) {
      brojac_C_redovi++;
      brojac_C_kolone = 0;
      }
  }
  
 // ofstream FILE("matrica.txt");
  cout<<endl<<endl;
   cout<<"M_M:"<<endl<<endl;
    for(int i = 0;i<M_M.size();i++) {
        
        for(int j = 0;j<M_M[i].size();j++) {
             
          
            cout<<M_M[i][j]<<" ";
            
        }
        
        cout<<endl;
        
    }
    
    
   
    
    
    
    
  /*vector<vector<double>> M_12 = { { 5, -2, 2, 7 },
                    { 1, 0, 0, 3 },
                    { -3, 1, 5, 0 },
                    { 3, -1, -9, 4 } };
                    
     */               
  vector<double> V(M_M[0].size(),0.0);
  vector<vector<double>> M_help(M_M.size(),V);
  
 
  vector<vector<double>> M_inv(M_M.size(),V);
  
  
 // double Matrica[N][N];
    adjungovana(M_M, M_help);

      
   // adjungovana(M_M,M_help);
  
   cout<<endl<<endl;
    cout<<"adjungovana:"<<endl<<endl;
    for(int i = 0;i<M_help.size();i++) {
        
        for(int j = 0;j<M_help[i].size();j++) {
            
            cout<<M_help[i][j]<<" ";
            
        }
        
        cout<<endl;
        
    }
    
        if (inverzna(M_M, M_inv)) {
            cout<<endl<<endl;
    cout<<"inverzna:"<<endl<<endl;
    for(int i = 0;i<M_inv.size();i++) {
        
        for(int j = 0;j<M_inv[i].size();j++) {
            
            cout<<M_inv[i][j]<<" ";
            
        }
        
        cout<<endl;
       
    }
        }
        
    
    vector<double> t_help(1,0);    
    vector<vector<double>> z_h(z.size(),t_help);
  //cout<<"zh: "<<endl<<endl;
    int counter_2 = 0;
    for(int i = 0;i<z_h.size();i++) {
        
        for(int j = 0;j<z_h[0].size();j++) {
            
            z_h[i][j] = z[counter_2] - h[counter_2];
       //     cout<<z_h[i][j];
        }
        cout<<endl;
        
        counter_2++;
    }
    
    vector<vector<double>> H_T_W_z_h = mnozenje_matrica(H_T_W,z_h);
     cout<<"H_T_W_z_h: "<<endl<<endl;
    for(int i = 0;i<H_T_W_z_h.size();i++) {
        
        for(int j = 0;j<H_T_W_z_h[i].size();j++) {
            
            cout<<H_T_W_z_h[i][j]<<" ";
        }
        cout<<endl;
    }
    
    vector<double> c1;
    vector<double> c2;
    double el_1 = (abs(Y[1][1])*abs(V2)*cos(arg(Y[1][1])+delta2) + abs(Y[1][2])*abs(V3)*cos(arg(Y[1][2])+delta3)+ abs(Y[1][3])*abs(V4)*cos(arg(Y[1][3]+delta4)));
    double el_2 = (abs(Y[1][1])*abs(V2)*sin(arg(Y[1][1])+delta2) + abs(Y[1][2])*abs(V3)*sin(arg(Y[1][2])+delta3)+ abs(Y[1][3])*abs(V4)*sin(arg(Y[1][3]+delta4)));
    
    c1.push_back(-1*el_1);
    c2.push_back(-1*el_2);
    
    
    
    H_T_W_z_h.push_back(c1);
    H_T_W_z_h.push_back(c2);
    
    
     cout<<"H_T_W_z_h: "<<endl<<endl;
    for(int i = 0;i<H_T_W_z_h.size();i++) {
        
        for(int j = 0;j<H_T_W_z_h[i].size();j++) {
            cout<<H_T_W_z_h[i][j]<<" ";
      //      FILE<<H_T_W_z_h[i][j]<<" ";
        }
        cout<<endl;
       
    }
    
    
    
    vector<vector<double>> rez_GLAVNI = mnozenje_matrica(M_inv,H_T_W_z_h);
    
    
     cout<<endl<<endl;
    cout<<"glavna:"<<endl<<endl;
    vector<double> konv;
    for(int i = 0;i<rez_GLAVNI.size();i++) {
        
           
        for(int j = 0;j<rez_GLAVNI[i].size();j++) {
            
           
        
           
            if(i == 0) {
               
                 konv.push_back(abs(delta2-rez_GLAVNI[i][j]));
                cout<<"delta2: "<<delta2<<endl;
                cout<<"delta2+1: "<<rez_GLAVNI[i][j]<<endl;
            }
            
            if(i == 1) {
                 
                 konv.push_back(abs(delta3-rez_GLAVNI[i][j]));
                  cout<<"delta3: "<<delta3<<endl;
                cout<<"delta3+1: "<<rez_GLAVNI[i][j]<<endl;
            }
            
             if(i == 2) {
                   
                 konv.push_back(abs(delta4-rez_GLAVNI[i][j]));
                  cout<<"delta4: "<<delta4<<endl;
                cout<<"delta4+1: "<<rez_GLAVNI[i][j]<<endl;
            }
             if(i == 3) {
                 
                 konv.push_back(abs(V1-rez_GLAVNI[i][j]));
                  cout<<"v1: "<<V1<<endl;
                cout<<"v1+1: "<<rez_GLAVNI[i][j]<<endl;
                 
            }
            
            if(i == 4) {
                 
                 konv.push_back(abs(V2-rez_GLAVNI[i][j]));
                 cout<<"v2: "<<V2<<endl;
                cout<<"v2+1: "<<rez_GLAVNI[i][j]<<endl;
            }
            
            if(i == 5) {
                 
                 konv.push_back(abs(V3-rez_GLAVNI[i][j]));
                 cout<<"v3: "<<V3<<endl;
                cout<<"v3+1: "<<rez_GLAVNI[i][j]<<endl;
            }
            
            if(i == 6) {
                 
                konv.push_back(abs(V4-rez_GLAVNI[i][j]));
                cout<<"v4: "<<V4<<endl;
                cout<<"v4+1: "<<rez_GLAVNI[i][j]<<endl;
            }
            
            
            
        }
        
        cout<<endl;
        
    }
    cout<<"konv: "<<endl<<endl;
    for(int i=0;i<7;i++) {
        cout<<konv[i]<<" ";
    }
    cout<<endl;
    
    
    if(counter_iteracija == 3) {
        
      
        rez_vek.push_back(delta2);
         rez_vek.push_back(delta3);
          rez_vek.push_back(delta4);
           rez_vek.push_back(V1);
            rez_vek.push_back(V2);
             rez_vek.push_back(V3);
              rez_vek.push_back(V4);
    }
    
   double deltax =  *max_element(konv.begin(), konv.end());
   cout<<"deltax: "<<deltax<<endl<<endl;
   
   if( deltax > epsilon) {
       cout<<counter_iteracija+1<<". iteracija "<<endl;
       delta2 += rez_GLAVNI[0][0];
       delta3 += rez_GLAVNI[1][0];
       delta4 += rez_GLAVNI[2][0];
       V1 += rez_GLAVNI[3][0];
       V2 += rez_GLAVNI[4][0];
       V3 += rez_GLAVNI[5][0];
       V4 += rez_GLAVNI[6][0];
       counter_iteracija++;
       continue;
       
   }
   else {
       cout<<endl<<endl;
       
       break;
   }
    
    
    }
    
    
    
  cout<<endl<<endl<<endl;
    
    
  cout<<"Program je završio nakon: "<<counter_iteracija-1<<". iteracije";
  
  cout<<endl<<"Rezultat: "<<endl;
  for(int i=0;i<7;i++) {
      cout<<rez_vek[i]<<" ";
  }

    
    return 0;
}








