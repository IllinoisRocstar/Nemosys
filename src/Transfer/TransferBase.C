#include <TransferBase.H>
#include <FETransfer.H>

TransferBase* TransferBase::Create(std::string method, meshBase* _source, meshBase* _target)
{
  if (!method.compare("Consistent Interpolation"))
  {
    FETransfer* transobj = new FETransfer( _source , _target);
    return transobj; 
  }
  else
  {
    std::cout << "Method " << method << " is not supported" << std::endl;
    std::cout << "Supported methods are: " << std::endl
              << "1) Consistent Interpolation" << std::endl;
    exit(1);
  }  
}


