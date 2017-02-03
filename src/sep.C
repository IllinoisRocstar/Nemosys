/// 
/// @file
/// @ingroup irproject_group
/// @brief Main for example serial program.
/// @author Mike Campbell (mtcampbe@illinois.edu)
/// @date 
/// 
#include "ExampleProgram.H"

typedef Nemosys::ExampleProgram::SEProgramType ProgramType;

int main(int argc,char *argv[])
{
  return(Nemosys::ExampleProgram::Driver<ProgramType>(argc,argv));
}
