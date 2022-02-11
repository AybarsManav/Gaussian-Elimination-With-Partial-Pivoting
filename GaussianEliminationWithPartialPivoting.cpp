#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
double** invertMartix(double** a);

//This function returns the one norm of an nxn matrix 
double oneNorm(double**a, int n)
{
    double abs_max_column_sum=0;
    for (int j = 0; j < n ; j++)
    {
        double sum = 0;
        for(int i= 0; i<n; i++) //Summing the absolute value of elements in column j
        {
            sum += abs(a[i][j]);
        }
        if (sum > abs_max_column_sum) //comparing it to abs_max_column_sum
        {
            abs_max_column_sum = sum;
        }
        
    }
    return abs_max_column_sum;
}

//This function returns the infinite norm of an nxn matrix
double infiniteNorm(double** a, int n)
{
    double abs_max_row_sum=0;
    for (int i = 0; i < n; i++)
    {
        double sum = 0;
        for (int j = 0; j < n; j++)
        {
            sum += abs(a[i][j]);
        }
        if (sum > abs_max_row_sum)
        {
            abs_max_row_sum = sum;
        } 
    }
    return abs_max_row_sum;
}

//This function returns a condition vector of condition numbers of 1 and infinitiy respectively.
double* findCondition(double** a)
{   
  
    double** inv_a = new double*[2];
    for (int i = 0; i < 2; i++)
    {
        inv_a[i] = new double[2];
    }
    
    inv_a = invertMartix(a);
    
    double* conditionNumber= new double[2];
    conditionNumber[0] = oneNorm(a,2)*oneNorm(inv_a,2);
    conditionNumber[1] = infiniteNorm(a,2)*infiniteNorm(inv_a,2);
    return conditionNumber;
}

//This function returns the determinant of the given 2x2 matrix.
double findDeterminant(double** a)
{
    return a[0][0] * a[1][1] - a[0][1] * a[1][0];
}

//This function inverts a 2x2 matrix and returns it.
double** invertMartix(double** a)
{
    //if we have a 2x2 matrix to invert it we calculate swap the positions of a,d
    //and multiply b,c with -1 then divide it into the original matrix's determinant.
    double determinant = findDeterminant(a);
    
    //Initializing the dynamic matrix inv_a;
    double **inv_a = new double*[2];
    for (int i = 0; i < 2; i++)
    {
        inv_a[i] = new double[2];
    }
    

    inv_a[0][0] = a[1][1]/determinant;
    inv_a[1][1] = a[0][0]/determinant;
    inv_a[0][1] = -a[0][1]/determinant;
    inv_a[1][0] = -a[0][1]/determinant;
    
    return inv_a;
}
//Gives us the index of the row of A with the greatest abs value starting from row k,
//We start from row k because every time we do an elimiation in a column, next time doing elimination
//we look at a submatrix, which means we do not look at the previous rows if they have larger elements since 
//they cannot be our pivots.
int findGreatestEntryIndex(double** a, int column, int n)
{  
    int max_index = column;
    for(int i=column ; i<n;i++)
    {
        if(abs(a[i][column]) > abs(a[max_index][column]))
        {
            max_index = i;
        }
    }
    return max_index;
}
//This function takes the adresses of the elements in matrix A and vector b as well as the
//rows to be interchanged also the dimension of the matrix.
void swapRows(double** a, double* b, int row1 , int row2,int n){
    double temp[n][n];
 
    double tempb;
    for (int i = 0; i < n; i++)
    {
        temp[row1][i] = a[row1][i];
        a[row1][i] = a[row2][i];
        a[row2][i] = temp[row1][i]; //We have swapped row1 and row2 in matrix A
    }
    tempb = b[row1];
    b[row1] = b[row2];
    b[row2] = tempb;       
}

int main(int argc, char** argv)
    {
    ifstream fileA(argv[1]);
    ifstream fileb(argv[2]); //We are creating our ifstream objects and linking them with files A.txt and b.txt
    ofstream filex("x.txt"); // We also create our ofstream object and link it with the file x.txt

    string line = ""; //We initialize a string called lines to see read the files.
    int n=0;          //The number of lines we read from the files will mean that we will be dealing with nxn matrices.

    while(getline(fileb,line)) // Counting the number of the lines in fileb will tell us how the dimensions of the matrix A.
    {
        n++;
    }
   

	fileb.clear(); //We clear all of the flags after reading the file
	fileb.seekg(0,ios::beg);//Because we have already read from the fileb we set the reading location to beginning of the file
   //for further reading it

    double *b = new double[n]; // We declare our vectos x and b dynamically.
    double *x = new double[n]; 

    for (int i = 0; i < n; i++)
    {
        fileb >> b[i]; //We read the ith element of the vector b from the file
        x[i] = 0;   //And set every element of x to zero.
    }

    double **A = new double*[n]; // We declare we will have an nxn matrix, here A[i][j] is a value in the matrix which is equivalent to *(*(A+i)+j) 
    for (int i = 0; i < n; i++) // We have now created the nxn matrix A
    {
        A[i] = new double[n];
    }
    
    for (int i = 0; i < n; i++)// Now we will initialize the matrix A with nested for loops
    {
        for (int j = 0; j < n; j++)
        {
            fileA >> A[i][j]; // For easiness, i = row i, j = column j
        }
        
    }
    
    fileb.close(); //Closing the input files.
	fileA.close();

    //We are obliged to compute the condition numbers of 2x2 matrices.
    // 1-norm is the absolute maximum column sum of the matrix
    // infinite-norm is the absolute minimum row sum of the matrix
    //cond(A) = Norm(A)*Norm(A^-1)
    
    //If the matrix is 2x2 we need to calculate the condition numbers and print it out to terminal.
    if(n==2)
    {
        double condition[2];
        for (int i = 0; i < 2; i++)
        {
            condition[i] = findCondition(A)[i];
        }
        cout << "The condition number of 1 is " << condition[0] << endl;
        cout << "The condition number of infinity is " << condition[1] << endl;

        if(isinf(condition[0]) || condition[0] > 1.0e+08){ //Checking if condition numbers are too high, we assume the matrix is singular.
            cout<<"The matrix is singular."<<endl;
            return 0; //We quit if the matrix is singular.
        }
    }

    //Gaussian elimination part
 	int max_index;

    for (int column = 0; column < n; column++) 
    {
        max_index = findGreatestEntryIndex(A, column, n); //We got the index of max element.

        if(max_index != column)
        {
            swapRows(A, b, column, max_index, n); //We swapped the row having the max element with (column)th row.
        }

        if(abs(A[column][column]) < pow(2,-10))//We assume EMech is 2^-10
        {
            cout <<"The matrix A is singular"<<endl; // If the pivot is smallar then EMech, we say the matrix is singular.
     
            return 0;
        }


        //Now we have our pivot in the specific column indicated by outer for. 
        //We must find the multiplicative factors for each entry below the pivot.
        double factor;
        for (int i = column+1; i < n ; i++)//We start from (column+1)th row to get the entries below.
        {
            factor = A[i][column]/A[column][column]; //we get the factor to make A[i][column] 0.
            for (int j = column+1; j < n; j++)
            {
                A[i][j] -= factor* A[column][j]; //Now we have subtracted factor* row number (column) from row i;
               
            }
            b[i] -= factor*b[column]; // We 
            A[i][column] = 0; //Since the floating point arithmetic is not always giving us 0 after such an operation
            //we know A[i][column] -= factor* A[column][column] gives a number very close to 0, lower than EMach
            //We assume all the entries below A[column][column] is 0.

        }
        
    }
  
     //Back substitution part
    for (int i = n-1; i >= 0; i--) //we start from row n-1
    {   

        for (int j = n-1; i < j; j--) //since we always subtract the rightest entries, j>i is always true.
        {
            b[i] = b[i]- A[i][j]*x[j]; // we subtract the known x values times the corresponding element of the matrix
        }

        x[i] = b[i] / A[i][i]; //after subtraction we divide to the matrix element
        
    }

    //Now we can write the result vector x to file x.txt
    for (int i = 0; i < n; i++)
    {
        filex << x[i] << "\n";
    }
    
    


	filex.close();
	
	
	return 0;
}
