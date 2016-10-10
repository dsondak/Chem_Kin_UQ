// Testing use of grvy for reading input from a file 
//
// Author: dsondak@ices.utexas.edu (David Sondak)

#include "grvy.h"
#include <iostream>
#include "read_input.h"

using namespace GRVY;

GRVY_Input_Class iparse;

// 
int ReadInt() 
{

    int int_in;

    if (! iparse.Open("../tests/input.txt") )
    {
       int_in = 0;
    }

    iparse.Read_Var("int_in", &int_in);

    iparse.Close();

    return int_in;
}

double ReadDouble()
{
    double doub_in;
    if (! iparse.Open("../tests/input.txt") )
    {
       doub_in = 0;
    }

    iparse.Read_Var("doub_in", &doub_in);

    iparse.Close();

    return doub_in;
}

std::string ReadString()
{
    std::string string_in;
    if (! iparse.Open("../tests/input.txt") )
    {
       string_in = "Failed";
    }

    iparse.Read_Var("string_in", &string_in);

    iparse.Close();

    return string_in;
}
