#include "TransferBase.H"

#include "FETransfer.H"

TransferBase *TransferBase::Create(const std::string &method,
                                   meshBase *_source,
                                   meshBase *_target)
{
  if (method == "Consistent Interpolation")
  {
    auto *transobj = new FETransfer(_source, _target);
    return transobj;
  }
  else
  {
    std::cerr << "Method " << method << " is not supported\n";
    std::cerr << "Supported methods are: \n"
              << "1) Consistent Interpolation" << std::endl;
    exit(1);
  }
}
