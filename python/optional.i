// Also see https://github.com/swig/swig/issues/1307

// PyOptType is a dummy name to tell SWIG about the template
// OptTemplate should be a C++ template that has one type parameter
// OptTemplate<Type> should have the same interface as std::optional<Type>
%define OPTIONAL_TEMPLATE(PyOptType, OptTemplate, Type)

%template(PyOptType) OptTemplate<Type>;

%typemap(in, implicitconv=1) OptTemplate<Type> {
  auto &temp = $1;
  if ($input != Py_None) {
    Type $1;
    $typemap(in, Type)
    temp = std::move($1);
  } else {
    temp = {};
  }
}

%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER) OptTemplate<Type> {
  if ($input == Py_None) {
    $1 = true;
  } else {
    $typemap(typecheck, Type)
  }
}

%typemap(in, implicitconv=1) const OptTemplate<Type> & (OptTemplate<Type> temp) {
  if ($input != Py_None) {
    Type $1;
    $typemap(in, Type)
    temp = std::move($1);
  } else {
    temp = {};
  }
  $1 = &temp;
}

%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER) const OptTemplate<Type> & {
  if ($input == Py_None) {
    $1 = true;
  } else {
    $typemap(typecheck, const Type &)
  }
}

%typemap(out) OptTemplate<Type> {
  if ($1.has_value()) {
    auto& temp = $1;
    {
      auto &$1 = temp.value();
      $typemap(out, Type)
    }
  } else {
    $result = Py_None;
    Py_INCREF(Py_None);
  }
}

%typemap(out) OptTemplate<Type> * {
  if ($1->has_value()) {
    auto& temp = $1;
    {
      auto $1 = &temp->value();
      $typemap(out, Type *)
    }
  } else {
    $result = Py_None;
    Py_INCREF(Py_None);
  }
}

%typemap(out) const OptTemplate<Type> & {
  if ($1->has_value()) {
    auto& temp = $1;
    {
      auto $1 = &temp->value();
      $typemap(out, const Type &)
    }
  } else {
    $result = Py_None;
    Py_INCREF(Py_None);
  }
}

%typemap(out) OptTemplate<Type> & {
  if ($1->has_value()) {
    auto& temp = $1;
    {
      auto $1 = &temp->value();
      $typemap(out, Type &)
    }
  } else {
    $result = Py_None;
    Py_INCREF(Py_None);
  }
}

%enddef

// For immutable types, naturalvar instructs SWIG to treat public data members of that type as if they are
// read/written by "const &" getters/setters
%define OPTIONAL_TEMPLATE_IMMUTABLE(PyOptType, OptTemplate, Type)
%naturalvar OptTemplate<Type>;
OPTIONAL_TEMPLATE(%arg(PyOptType), %arg(OptTemplate), %arg(Type));
%enddef
