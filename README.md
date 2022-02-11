# Gaussian-Elimination-With-Partial-Pivoting
 
 In this project, I implemented the Gaussian elimination algortihm with partial pivoting togetger with backwards substitution to solve system like Ax=b
 where A is an nxn square matrix.
 
 The program has two command line arguments for the parameters:
 - The first argument is the name of the file we read the A matrix from.
 - The second argument is the name of the file we read the b vector from.
 
 Each line represents a row.

The program prints out an error message and quits if it detecs that A is singular, in this code the singularity condition for general nxn matrices is decided so that if a pivot is smallar than a predetermined number, the matrix is considered singular. On the other hand, the program also prints out the condition numbers of 2x2 matrices, consequently another singularity metric for such matrices is having high condition numbers.

The content of the x vector is saved in a txt file named "x.txt"

