#include <TransferBase.H>
#include <FETransfer.H>

TransferBase* TransferBase::Create(std::string method, meshBase* _source, meshBase* _target)
{
  if (!method.compare("FE"))
  {
    FETransfer* transobj = new FETransfer( _source , _target);
    return transobj; 
  }  
}


