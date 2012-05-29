/*
 *  GeneralFunctions.cpp
 *  Segments
 *
 *  Created by Michel Koskas on 22/08/11.
 *  Copyright 2011 INRA, INA. All rights reserved.
 *
 */

#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include "GeneralFunctions.h"
#include "MyVector.h"


// This function says whether the character x is a number or not
bool IsDigit(char &x)
{
  if (x < '0')
    return false;
  if (x > '9')
    return false;
  return true;
}

bool ToNext(char *Buffer, int &BuffIndex, int BufferSize, char Separator, char Terminator)
{
  while ((BuffIndex < BufferSize) && (Buffer[BuffIndex] != Separator) && (Buffer[BuffIndex] != Terminator))
    BuffIndex++;
  if (BuffIndex == BufferSize)
  {
    std::cerr << "Can NOT drop this element (Buffer finished too soon.)" << std::endl;
    std::cerr << "Exiting with errcode 121." << std::endl;
    exit(121);
  }
  BuffIndex++;
	return true;
}


int GetRandomNumber(int MinValue, int MaxValue)
{
  long int x = rand();
  double y = (MaxValue - MinValue + 1) * ((double) x) / RAND_MAX;
  int Res = MinValue + y;
  return Res;
}

MyVector<int> IntersectLists(const MyVector<int> &A, const MyVector<int> &B)
{
	MyVector<int> Res;
	MyVector<int>::const_iterator IA = A.begin(), IB = B.begin();
	while ((IA != A.end()) && (IB != B.end()))
		if (*IA < *IB)
			IA++;
		else if (*IB < *IA)
			IB++;
		else
		{
			Res.push_back(*IA);
			IA++;
			IB++;
		}
	return Res;
}

MyVector<int> GetBreakpoints(int k, int n, int** M)
{
  MyVector<int> TheBreakpoints;
  if (k>1)
    {
      int Prec = M[k-1][n-1];
      TheBreakpoints.push_back(Prec+1);
      if (k>2)
	for (int i=(k-2); i>=1; i--)
	  {
	    TheBreakpoints.push_back(M[i][Prec]+1);
	    Prec = M[i][Prec];
	  }
    }
  TheBreakpoints.push_back(0);
  TheBreakpoints.reverse();
  TheBreakpoints.push_back(n);
  TheBreakpoints.sort();
  return TheBreakpoints;
}

MyVector<double> GetParameters(int k, int n, int** M, double** Par)
{
  MyVector<double> Parameters;
  Parameters.push_back(Par[k-1][n-1]);
  if (k>1)
    {
      int Prec = M[k-1][n-1];
      if (k>2)
	    for (int i=(k-2); i>=1; i--)
	      {
	        Parameters.push_back(Par[i][Prec]);
	        Prec = M[i][Prec];
	      }
      Parameters.push_back(Par[0][Prec]);
    }
  Parameters.reverse();
  return Parameters;
}

void WriteAllWithTime(int N, int K, int nbODPA, int nbrupidentik, double Relative, double tPDPA, double tCart, const char* FileName)
{
  std::ofstream MyResults;
  MyResults.open (FileName, std::ios_base::app);
  if (!MyResults.is_open())
  {
	std::cerr << "I Can NOT open the file " << MyResults << " and I'm leaving this mess." << std::endl;
	exit(1465);
  }
  MyResults.seekp(0,std::ios_base::end);
  MyResults<< N << '\t' << K << '\t'<< nbODPA << '\t'<< nbrupidentik << '\t' << Relative << '\t' << tPDPA << '\t' << tCart  << '\t' << '\n' ;
  MyResults.close();
}

void WriteTypeWithTime(int choice, int N, int K, int nbODPA, double tPDPA, double tODPA, const char* FileName)
{
  std::ofstream MyResults;
  MyResults.open (FileName, std::ios_base::app);
  if (!MyResults.is_open())
  {
	std::cerr << "I Can NOT open the file " << MyResults << " and I'm leaving this mess." << std::endl;
	exit(1465);
  }
  MyResults.seekp(0,std::ios_base::end);
  MyResults << choice << '\t' << N << '\t' << K << '\t' << nbODPA << '\t'  << tPDPA << '\t' << tODPA << '\n' ;
  MyResults.close();
}

