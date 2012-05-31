/*
 *  GeneralFunctions.h
 *  Segments
 *
 *  Created by Michel Koskas on 22/08/11.
 *  Copyright 2011 INRA, INA. All rights reserved.
 *
 */
#ifndef GeneralFunctions_h_
#define GeneralFunctions_h_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "Constants.h"
#include "MyVector.h"


bool IsDigit(char &x);
bool ToNext(char *Buffer, int &BuffIndex, int BufferSize, char Separator=Separateur, char Terminator=FinDeLigne);
int GetRandomNumber(int MinValue, int MaxValue);


template <typename T>
void Swap(T &A, T &B)
{
	T C = A;
	A = B;
	B = C;
}

// In this function, Buffer is a table of characters, BuffIndex is the current reading location and BufferSize its global size.
// The next integer is read and its value is stocked in Res.

template <typename DataTypeName>
bool NextNumber(char *Buffer, int &BuffIndex, int &BufferSize, DataTypeName &Res, bool &VuRetourChariot, char Separator=Separateur, char Terminator=FinDeLigne)
{
  Res = 0;
  VuRetourChariot = false;
  while ((BuffIndex < BufferSize) && (!IsDigit(Buffer[BuffIndex]) && (Buffer[BuffIndex] != '-')))
  {
    if (Buffer[BuffIndex] == '\n')
      VuRetourChariot = true;
    BuffIndex++;
  }
  if (BuffIndex == BufferSize)
  {
    Res = -1;
    return false;
  }
	int Signe = 1;
	if (Buffer[BuffIndex] == '-')
	{
		Signe = -1;
		BuffIndex++;
	}
  Res = 0;
  while ((BuffIndex < BufferSize) && (IsDigit(Buffer[BuffIndex])))
    Res = 10 * Res + Buffer[BuffIndex++] - '0';
  if (Buffer[BuffIndex] == '.')
  {
    BuffIndex++;
    double P = 1;
    while ((BuffIndex < BufferSize) && (IsDigit(Buffer[BuffIndex])))
    {
      P /= 10;
      Res += (Buffer[BuffIndex++] - '0') * P;
    }
  }
	Res *= Signe;
	int FormerBufferIndex = BuffIndex;
	ToNext(Buffer, BuffIndex, BufferSize, Separator, Terminator);
	for (int i = FormerBufferIndex; i < BuffIndex; i++)
		if (Buffer[i] == Terminator)
		{
			VuRetourChariot = true;
			break;
		}
  return true;
}

#endif
