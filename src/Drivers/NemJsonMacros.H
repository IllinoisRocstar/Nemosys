#ifndef NEMOSYS_NEMJSONMACROS_H_
#define NEMOSYS_NEMJSONMACROS_H_
#include <jsoncons/json_traits_macros.hpp>

#define NEM_JSON_UNPACK(...) __VA_ARGS__

// Extending "Mode": NEM_JSON_RDONLY_OVERRIDE and NEM_JSON_RDWR_OVERRIDE act
// like JSONCONS_RDONLY and JSONCONS_RDWR in the "as" methods, but are no-ops in
// the "to_json" methods, as it is assumed that the parent to_json will
// serialize the member. Take care to use the same serialized name!
#define NEM_JSON_RDONLY_OVERRIDE(X)
#define NEM_JSON_RDWR_OVERRIDE(X) X
#define JSONCONS_RDONLY_TO_JSON_MODE(X) X
#define JSONCONS_RDWR_TO_JSON_MODE(X) X
#define NEM_JSON_RDONLY_OVERRIDE_TO_JSON_MODE(X)
#define NEM_JSON_RDWR_OVERRIDE_TO_JSON_MODE(X)

#define NEM_JSON_MEMBER_NAME_IS(P1, P2, P3, Seq, Count) NEM_JSON_MEMBER_NAME_IS_LAST(P1, P2, P3, Seq, Count)
#define NEM_JSON_MEMBER_NAME_IS_LAST(P1, P2, P3, Seq, Count) if ((num_params-Count) < num_mandatory_params1 && JSONCONS_EXPAND(JSONCONS_CONCAT(NEM_JSON_MEMBER_NAME_IS_,JSONCONS_NARGS Seq) Seq)
#define NEM_JSON_MEMBER_NAME_IS_2(Member, Name) !ajson.contains(Name)) return false;
#define NEM_JSON_MEMBER_NAME_IS_3(Member, Name, Mode) !ajson.contains(Name)) return false; \
    if (ajson.contains(Name) && !ajson.at(Name).template is<typename std::decay<decltype(((value_type*)nullptr)->Member)>::type>()) return false;
#define NEM_JSON_MEMBER_NAME_IS_4(Member, Name, Mode, Match) NEM_JSON_MEMBER_NAME_IS_6(Member, Name, Mode, Match, , )
#define NEM_JSON_MEMBER_NAME_IS_5(Member, Name, Mode, Match, Into) NEM_JSON_MEMBER_NAME_IS_6(Member, Name, Mode, Match, Into, )
#define NEM_JSON_MEMBER_NAME_IS_6(Member, Name, Mode, Match, Into, From) !ajson.contains(Name)) return false; \
    JSONCONS_TRY{if (ajson.contains(Name) && !Match(ajson.at(Name).template as<typename std::decay<decltype(Into(((value_type*)nullptr)->Member))>::type>())) return false;} \
    JSONCONS_CATCH(...) {return false;}

#define NEM_JSON_N_MEMBER_NAME_AS(P1, P2, P3, Seq, Count) NEM_JSON_N_MEMBER_NAME_AS_LAST(P1, P2, P3, Seq, Count)
#define NEM_JSON_N_MEMBER_NAME_AS_LAST(P1, P2, P3, Seq, Count) JSONCONS_EXPAND(JSONCONS_CONCAT(NEM_JSON_N_MEMBER_NAME_AS_,JSONCONS_NARGS Seq) Seq)
#define NEM_JSON_N_MEMBER_NAME_AS_2(Member, Name) \
    if (ajson.contains(Name)) {json_traits_helper<Json>::set_udt_member(ajson,Name,aval->Member);}
#define NEM_JSON_N_MEMBER_NAME_AS_3(Member, Name, Mode) Mode(NEM_JSON_N_MEMBER_NAME_AS_2(Member, Name))
#define NEM_JSON_N_MEMBER_NAME_AS_4(Member, Name, Mode, Match) \
    Mode(if (ajson.contains(Name)) {json_traits_helper<Json>::set_udt_member(ajson,Name,aval->Member);})
#define NEM_JSON_N_MEMBER_NAME_AS_5(Member, Name, Mode, Match, Into) \
    Mode(if (ajson.contains(Name)) {json_traits_helper<Json>::template set_udt_member<typename std::decay<decltype(Into(((value_type*)nullptr)->Member))>::type>(ajson,Name,aval->Member);})
#define NEM_JSON_N_MEMBER_NAME_AS_6(Member, Name, Mode, Match, Into, From) \
    Mode(if (ajson.contains(Name)) {json_traits_helper<Json>::template set_udt_member<typename std::decay<decltype(Into(((value_type*)nullptr)->Member))>::type>(ajson,Name,From,aval->Member);})

#define NEM_JSON_N_MEMBER_NAME_TO_JSON(P1, P2, P3, Seq, Count) NEM_JSON_N_MEMBER_NAME_TO_JSON_LAST(P1, P2, P3, Seq, Count)
#define NEM_JSON_N_MEMBER_NAME_TO_JSON_LAST(P1, P2, P3, Seq, Count) if ((num_params-Count) < num_mandatory_params2) JSONCONS_EXPAND(JSONCONS_CONCAT(NEM_JSON_N_MEMBER_NAME_TO_JSON_,JSONCONS_NARGS Seq) Seq)
#define NEM_JSON_N_MEMBER_NAME_TO_JSON_2(Member, Name) \
  {ajson.try_emplace(Name, aval->Member);} \
else \
  {json_traits_helper<Json>::set_optional_json_member(Name, aval->Member, ajson);}
#define NEM_JSON_N_MEMBER_NAME_TO_JSON_3(Member, Name, Mode) NEM_JSON_N_MEMBER_NAME_TO_JSON_2(Member, Name)
#define NEM_JSON_N_MEMBER_NAME_TO_JSON_4(Member, Name, Mode, Match) NEM_JSON_N_MEMBER_NAME_TO_JSON_6(Member, Name, Mode, Match,,)
#define NEM_JSON_N_MEMBER_NAME_TO_JSON_5(Member, Name, Mode, Match, Into) NEM_JSON_N_MEMBER_NAME_TO_JSON_6(Member, Name, Mode, Match, Into, )
#define NEM_JSON_N_MEMBER_NAME_TO_JSON_6(Member, Name, Mode, Match, Into, From) \
  {Mode ## _TO_JSON_MODE(ajson.try_emplace(Name, Into(aval->Member));)} \
else \
  {Mode ## _TO_JSON_MODE(json_traits_helper<Json>::set_optional_json_member(Name, Into(aval->Member), ajson);)}

#define NEM_JSON_N_MEMBER_NAME_TRAITS_FINAL(ValueType, ParentType, NumMandatoryParams, ...)  \
    namespace jsoncons \
{ \
    template<typename Json> \
    struct json_type_traits<Json, ValueType *> \
    { \
        using value_type = ValueType; \
        using allocator_type = typename Json::allocator_type; \
        using char_type = typename Json::char_type; \
        using string_view_type = typename Json::string_view_type; \
        constexpr static size_t num_params = JSONCONS_NARGS(__VA_ARGS__); \
        constexpr static size_t num_mandatory_params1 = NumMandatoryParams; \
        constexpr static size_t num_mandatory_params2 = NumMandatoryParams; \
        static bool is(const Json& ajson) noexcept \
        { \
            if (!ajson.template is<ParentType *>()) return false; \
            if (!ajson.is_object()) return false; \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_MEMBER_NAME_IS,,,, __VA_ARGS__)\
            return true; \
        } \
        static value_type *as(const Json& ajson) \
        { \
            if (!is(ajson)) JSONCONS_THROW(convert_error(convert_errc::conversion_failed, "Not a " # ValueType)); \
            value_type *aval = new value_type; \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_N_MEMBER_NAME_AS,,,, __VA_ARGS__) \
            json_type_traits<Json, ParentType *>::as_from_child(ajson, aval); \
            return aval; \
        } \
        static Json to_json(const value_type *const &aval, allocator_type alloc=allocator_type()) \
        { \
            Json ajson(json_object_arg, semantic_tag::none, alloc); \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_N_MEMBER_NAME_TO_JSON,,,, __VA_ARGS__) \
            json_type_traits<Json, ParentType *>::to_json_from_child(aval, ajson); \
            return ajson; \
        } \
    }; \
    template <> struct is_json_type_traits_declared<ValueType *> : public std::true_type {};\
} \
 /**/

#define NEM_JSON_POLYMORPHIC_AS_RAW_PTR(BaseClass, P2, P3, DerivedClass, Count) NEM_JSON_POLYMORPHIC_AS_RAW_PTR_LAST(BaseClass, P2, P3, DerivedClass, Count)
#define NEM_JSON_POLYMORPHIC_AS_RAW_PTR_LAST(BaseClass, P2, P3, DerivedClass, Count)  if (ajson.template is<DerivedClass *>()) return ajson.template as<DerivedClass *>();

#define NEM_JSON_POLYMORPHIC_TO_JSON(BaseClass, P2, P3, DerivedClass, Count) NEM_JSON_POLYMORPHIC_TO_JSON_LAST(BaseClass, P2, P3, DerivedClass, Count)
#define NEM_JSON_POLYMORPHIC_TO_JSON_LAST(BaseClass, P2, P3, DerivedClass, Count) { auto p = dynamic_cast<const DerivedClass *>(ptr); if (p) { return Json(const_cast<DerivedClass *>(p)); } }

#define NEM_JSON_N_MEMBER_NAME_TRAITS_BASE(BaseClass, ChildClasses, NumMandatoryParams, ...) \
namespace jsoncons { \
    template<class Json> \
    struct json_type_traits<Json, BaseClass *> { \
        using value_type = BaseClass; \
        using allocator_type = typename Json::allocator_type; \
        using char_type = typename Json::char_type; \
        using string_view_type = typename Json::string_view_type; \
        constexpr static size_t num_params = JSONCONS_NARGS(__VA_ARGS__); \
        constexpr static size_t num_mandatory_params1 = NumMandatoryParams; \
        constexpr static size_t num_mandatory_params2 = NumMandatoryParams; \
        static bool is(const Json& ajson) noexcept { \
            if (!ajson.is_object()) return false; \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_MEMBER_NAME_IS,,,, __VA_ARGS__)\
            return true; \
        } \
\
        static BaseClass *as(const Json& ajson) { \
            if (!ajson.is_object()) return nullptr; \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_POLYMORPHIC_AS_RAW_PTR, BaseClass,,, NEM_JSON_UNPACK ChildClasses) \
            JSONCONS_THROW(convert_error(convert_errc::conversion_failed, "Not a " # BaseClass)); \
        } \
\
        static void as_from_child(const Json &ajson, value_type * const&aval) { \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_N_MEMBER_NAME_AS,,,, __VA_ARGS__) \
        } \
\
        static Json to_json(const BaseClass * const& ptr) { \
            if (!ptr) {return Json::null();} \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_POLYMORPHIC_TO_JSON, BaseClass,,, NEM_JSON_UNPACK ChildClasses) \
            return Json::null(); \
        } \
\
        static void to_json_from_child(const BaseClass * const& aval, Json &ajson) { \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_N_MEMBER_NAME_TO_JSON,,,, __VA_ARGS__) \
        } \
    }; \
    template <> struct is_json_type_traits_declared<BaseClass *> : public std::true_type {}; \
}  \
  /**/

#define NEM_JSON_N_MEMBER_NAME_TRAITS_INTERMEDIATE(BaseClass, ChildClasses, ParentClass, NumMandatoryParams, ...) \
namespace jsoncons { \
    template<class Json> \
    struct json_type_traits<Json, BaseClass *> { \
        using value_type = BaseClass; \
        using allocator_type = typename Json::allocator_type; \
        using char_type = typename Json::char_type; \
        using string_view_type = typename Json::string_view_type; \
        constexpr static size_t num_params = JSONCONS_NARGS(__VA_ARGS__); \
        constexpr static size_t num_mandatory_params1 = NumMandatoryParams; \
        constexpr static size_t num_mandatory_params2 = NumMandatoryParams; \
        static bool is(const Json& ajson) noexcept { \
            if (!ajson.template is<ParentClass *>()) return false; \
            if (!ajson.is_object()) return false; \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_MEMBER_NAME_IS,,,, __VA_ARGS__)\
            return true; \
        } \
\
        static BaseClass * as(const Json& ajson) { \
            if (!ajson.is_object()) return nullptr; \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_POLYMORPHIC_AS_RAW_PTR, BaseClass,,, NEM_JSON_UNPACK ChildClasses) \
            JSONCONS_THROW(convert_error(convert_errc::conversion_failed, "Not a " # BaseClass)); \
        } \
\
        static void as_from_child(const Json &ajson, value_type * const &aval) { \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_N_MEMBER_NAME_AS,,,, __VA_ARGS__) \
            json_type_traits<Json, ParentClass *>::as_from_child(ajson, aval); \
        } \
\
        static Json to_json(const BaseClass * const& ptr) { \
            if (!ptr) {return Json::null();} \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_POLYMORPHIC_TO_JSON, BaseClass,,, NEM_JSON_UNPACK ChildClasses) \
            return Json::null(); \
        } \
\
        static void to_json_from_child(const BaseClass * const& aval, Json &ajson) { \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_N_MEMBER_NAME_TO_JSON,,,, __VA_ARGS__) \
            json_type_traits<Json, ParentClass *>::to_json_from_child(aval, ajson); \
        } \
    }; \
    template <> struct is_json_type_traits_declared<BaseClass *> : public std::true_type {}; \
}  \
  /**/

#define NEM_JSON_WRAP_SMART_PTR(BaseClass, PointerTempl) \
namespace jsoncons { \
    template<class Json> \
    struct json_type_traits<Json, PointerTempl<BaseClass>> { \
        using allocator_type = typename Json::allocator_type; \
\
        static bool is(const Json &ajson) noexcept { return ajson.template is<BaseClass *>(); }; \
\
        static PointerTempl<BaseClass> as(const Json &ajson) { \
            return PointerTempl<BaseClass>(ajson.template as<BaseClass *>()); \
        }; \
\
        static Json to_json(const PointerTempl<BaseClass> &aval) { \
            return Json(aval.get()); \
        } \
    }; \
\
    template<> \
    struct is_json_type_traits_declared<PointerTempl<BaseClass>> : public std::true_type { \
    }; \
} \
  /**/

#define NEM_JSON_SMART_PTR_VAL(BaseClass, PointerTempl) \
namespace jsoncons { \
    template<class Json> \
    struct json_type_traits<Json, PointerTempl<BaseClass>> { \
        using allocator_type = typename Json::allocator_type; \
\
        static bool is(const Json &ajson) noexcept { return ajson.template is<BaseClass>(); }; \
\
        static PointerTempl<BaseClass> as(const Json &ajson) { \
            return PointerTempl<BaseClass>(new BaseClass(ajson.template as<BaseClass>())); \
        }; \
\
        static Json to_json(const PointerTempl<BaseClass> &aval) { \
            return aval ? Json(*aval) : Json::null(); \
        } \
    }; \
\
    template<> \
    struct is_json_type_traits_declared<PointerTempl<BaseClass>> : public std::true_type { \
    }; \
} \
  /**/

#define NEM_JSON_GETTER_SETTER_NAME_IS(P1, P2, P3, Seq, Count) NEM_JSON_GETTER_SETTER_NAME_IS_LAST(P1, P2, P3, Seq, Count)
#define NEM_JSON_GETTER_SETTER_NAME_IS_LAST(P1, P2, P3, Seq, Count) if ((num_params-Count) < num_mandatory_params1 && JSONCONS_EXPAND(JSONCONS_CONCAT(NEM_JSON_GETTER_SETTER_NAME_IS_,JSONCONS_NARGS Seq) Seq)
#define NEM_JSON_GETTER_SETTER_NAME_IS_3(Getter, Setter, Name) !ajson.contains(Name)) return false;
#define NEM_JSON_GETTER_SETTER_NAME_IS_4(Getter, Setter, Name, Mode) !ajson.contains(Name)) return false; \
    if (ajson.contains(Name) && !ajson.at(Name).template is<typename std::decay<decltype(((value_type*)nullptr)->Getter())>::type>()) return false;
#define NEM_JSON_GETTER_SETTER_NAME_IS_5(Getter, Setter, Name, Mode, Match) NEM_JSON_GETTER_SETTER_NAME_IS_7(Getter, Setter, Name, Mode, Match,, )
#define NEM_JSON_GETTER_SETTER_NAME_IS_6(Getter, Setter, Name, Mode, Match, Into) NEM_JSON_GETTER_SETTER_NAME_IS_7(Getter, Setter, Name, Mode, Match, Into, )
#define NEM_JSON_GETTER_SETTER_NAME_IS_7(Getter, Setter, Name, Mode, Match, Into, From) !ajson.contains(Name)) return false; \
    JSONCONS_TRY{if (ajson.contains(Name) && !Match(ajson.at(Name).template as<typename std::decay<decltype(Into(((value_type*)nullptr)->Getter()))>::type>())) return false;} \
    JSONCONS_CATCH(...) {return false;}

#define NEM_JSON_N_GETTER_SETTER_NAME_AS(P1, P2, P3, Seq, Count) NEM_JSON_N_GETTER_SETTER_NAME_AS_LAST(P1, P2, P3, Seq, Count)
#define NEM_JSON_N_GETTER_SETTER_NAME_AS_LAST(P1, P2, P3, Seq, Count) JSONCONS_EXPAND(JSONCONS_CONCAT(NEM_JSON_N_GETTER_SETTER_NAME_AS_,JSONCONS_NARGS Seq) Seq)
#define NEM_JSON_N_GETTER_SETTER_NAME_AS_3(Getter, Setter, Name) if (ajson.contains(Name)) aval->Setter(ajson.at(Name).template as<typename std::decay<decltype(aval->Getter())>::type>());
#define NEM_JSON_N_GETTER_SETTER_NAME_AS_4(Getter, Setter, Name, Mode) Mode(NEM_JSON_N_GETTER_SETTER_NAME_AS_3(Getter, Setter, Name))
#define NEM_JSON_N_GETTER_SETTER_NAME_AS_5(Getter, Setter, Name, Mode, Match) NEM_JSON_N_GETTER_SETTER_NAME_AS_7(Getter, Setter, Name, Mode, Match, , )
#define NEM_JSON_N_GETTER_SETTER_NAME_AS_6(Getter, Setter, Name, Mode, Match, Into) NEM_JSON_N_GETTER_SETTER_NAME_AS_7(Getter, Setter, Name, Mode, Match, Into, )
#define NEM_JSON_N_GETTER_SETTER_NAME_AS_7(Getter, Setter, Name, Mode, Match, Into, From) Mode(if (ajson.contains(Name)) aval->Setter(From(ajson.at(Name).template as<typename std::decay<decltype(Into(aval->Getter()))>::type>()));)

#define NEM_JSON_N_GETTER_SETTER_NAME_TO_JSON(P1, P2, P3, Seq, Count) NEM_JSON_N_GETTER_SETTER_NAME_TO_JSON_LAST(P1, P2, P3, Seq, Count)
#define NEM_JSON_N_GETTER_SETTER_NAME_TO_JSON_LAST(P1, P2, P3, Seq, Count) JSONCONS_EXPAND(JSONCONS_CONCAT(NEM_JSON_N_GETTER_SETTER_NAME_TO_JSON_,JSONCONS_NARGS Seq) Seq)
#define NEM_JSON_N_GETTER_SETTER_NAME_TO_JSON_3(Getter, Setter, Name) ajson.try_emplace(Name, aval->Getter() );
#define NEM_JSON_N_GETTER_SETTER_NAME_TO_JSON_4(Getter, Setter, Name, Mode) NEM_JSON_N_GETTER_SETTER_NAME_TO_JSON_7(Getter, Setter, Name, Mode, , , )
#define NEM_JSON_N_GETTER_SETTER_NAME_TO_JSON_5(Getter, Setter, Name, Mode, Match) NEM_JSON_N_GETTER_SETTER_NAME_TO_JSON_7(Getter, Setter, Name, Mode, Match, , )
#define NEM_JSON_N_GETTER_SETTER_NAME_TO_JSON_6(Getter, Setter, Name, Mode, Match, Into) NEM_JSON_N_GETTER_SETTER_NAME_TO_JSON_7(Getter, Setter, Name, Mode, Match, Into, )
#define NEM_JSON_N_GETTER_SETTER_NAME_TO_JSON_7(Getter, Setter, Name, Mode, Match, Into, From) Mode ## _TO_JSON_MODE(ajson.try_emplace(Name, Into(aval->Getter()) );)

#define NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_BASE(BaseClass,ChildClasses,NumMandatoryParams, ...)  \
namespace jsoncons \
{ \
    template<typename Json> \
    struct json_type_traits<Json, BaseClass *> \
    { \
        using value_type = BaseClass; \
        using allocator_type = typename Json::allocator_type; \
        using char_type = typename Json::char_type; \
        using string_view_type = typename Json::string_view_type; \
        constexpr static size_t num_params = JSONCONS_NARGS(__VA_ARGS__); \
        constexpr static size_t num_mandatory_params1 = NumMandatoryParams; \
        constexpr static size_t num_mandatory_params2 = NumMandatoryParams; \
        static bool is(const Json& ajson) noexcept \
        { \
            if (!ajson.is_object()) return false; \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_GETTER_SETTER_NAME_IS,,,, __VA_ARGS__)\
            return true; \
        } \
        static value_type *as(const Json& ajson) { \
            if (!ajson.is_object()) return nullptr; \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_POLYMORPHIC_AS_RAW_PTR, BaseClass,,, NEM_JSON_UNPACK ChildClasses) \
            JSONCONS_THROW(convert_error(convert_errc::conversion_failed, "Not a " # BaseClass)); \
        } \
        static void as_from_child(const Json &ajson, value_type * const&aval) { \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_N_GETTER_SETTER_NAME_AS,,,, __VA_ARGS__) \
        } \
        static Json to_json(const BaseClass * const& ptr) { \
            if (!ptr) {return Json::null();} \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_POLYMORPHIC_TO_JSON, BaseClass,,, NEM_JSON_UNPACK ChildClasses) \
            return Json::null(); \
        } \
        static void to_json_from_child(const BaseClass * const& aval, Json &ajson) { \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_N_GETTER_SETTER_NAME_TO_JSON,,,, __VA_ARGS__) \
        } \
    }; \
    template <> struct is_json_type_traits_declared<BaseClass *> : public std::true_type {};\
} \
  /**/
#define NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_INTERMEDIATE(BaseClass, ChildClasses, ParentClass, NumMandatoryParams, ...) \
namespace jsoncons { \
    template<class Json> \
    struct json_type_traits<Json, BaseClass *> { \
        using value_type = BaseClass; \
        using allocator_type = typename Json::allocator_type; \
        using char_type = typename Json::char_type; \
        using string_view_type = typename Json::string_view_type; \
        constexpr static size_t num_params = JSONCONS_NARGS(__VA_ARGS__); \
        constexpr static size_t num_mandatory_params1 = NumMandatoryParams; \
        constexpr static size_t num_mandatory_params2 = NumMandatoryParams; \
        static bool is(const Json& ajson) noexcept { \
            if (!ajson.template is<ParentClass *>()) return false; \
            if (!ajson.is_object()) return false; \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_GETTER_SETTER_NAME_IS,,,, __VA_ARGS__)\
            return true; \
        } \
\
        static BaseClass * as(const Json& ajson) { \
            if (!ajson.is_object()) return nullptr; \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_POLYMORPHIC_AS_RAW_PTR, BaseClass,,, NEM_JSON_UNPACK ChildClasses) \
            JSONCONS_THROW(convert_error(convert_errc::conversion_failed, "Not a " # BaseClass)); \
        } \
\
        static void as_from_child(const Json &ajson, value_type * const &aval) { \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_N_GETTER_SETTER_NAME_AS,,,, __VA_ARGS__) \
            json_type_traits<Json, ParentClass *>::as_from_child(ajson, aval); \
        } \
\
        static Json to_json(const BaseClass * const& ptr) { \
            if (!ptr) {return Json::null();} \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_POLYMORPHIC_TO_JSON, BaseClass,,, NEM_JSON_UNPACK ChildClasses) \
            return Json::null(); \
        } \
\
        static void to_json_from_child(const BaseClass * const& aval, Json &ajson) { \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_N_GETTER_SETTER_NAME_TO_JSON,,,, __VA_ARGS__) \
            json_type_traits<Json, ParentClass *>::to_json_from_child(aval, ajson); \
        } \
    }; \
    template <> struct is_json_type_traits_declared<BaseClass *> : public std::true_type {}; \
}  \
  /**/

#define NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(ValueType,ParentType,NumMandatoryParams, ...)  \
namespace jsoncons \
{ \
    template<typename Json> \
    struct json_type_traits<Json, ValueType *> \
    { \
        using value_type = ValueType; \
        using allocator_type = typename Json::allocator_type; \
        using char_type = typename Json::char_type; \
        using string_view_type = typename Json::string_view_type; \
        constexpr static size_t num_params = JSONCONS_NARGS(__VA_ARGS__); \
        constexpr static size_t num_mandatory_params1 = NumMandatoryParams; \
        constexpr static size_t num_mandatory_params2 = NumMandatoryParams; \
        static bool is(const Json& ajson) noexcept \
        { \
            if (!ajson.template is<ParentType *>()) return false; \
            if (!ajson.is_object()) return false; \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_GETTER_SETTER_NAME_IS,,,, __VA_ARGS__)\
            return true; \
        } \
        static value_type *as(const Json& ajson) \
        { \
            if (!is(ajson)) JSONCONS_THROW(convert_error(convert_errc::conversion_failed, "Not a " # ValueType)); \
            value_type *aval = new value_type; \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_N_GETTER_SETTER_NAME_AS,,,, __VA_ARGS__) \
            json_type_traits<Json, ParentType *>::as_from_child(ajson, aval); \
            return aval; \
        } \
        static Json to_json(const value_type *const &aval, allocator_type alloc=allocator_type()) \
        { \
            Json ajson(json_object_arg, semantic_tag::none, alloc); \
            JSONCONS_VARIADIC_REP_N(NEM_JSON_N_GETTER_SETTER_NAME_TO_JSON,,,, __VA_ARGS__) \
            json_type_traits<Json, ParentType *>::to_json_from_child(aval, ajson); \
            return ajson; \
        } \
    }; \
    template <> struct is_json_type_traits_declared<ValueType *> : public std::true_type {};\
} \
  /**/

#define NEM_JSON_N_MEMBER_NAME_TO_JSON_VAL(P1, P2, P3, Seq, Count) NEM_JSON_N_MEMBER_NAME_TO_JSON_VAL_LAST(P1, P2, P3, Seq, Count)
#define NEM_JSON_N_MEMBER_NAME_TO_JSON_VAL_LAST(P1, P2, P3, Seq, Count) if ((num_params-Count) < num_mandatory_params2) JSONCONS_EXPAND(JSONCONS_CONCAT(NEM_JSON_N_MEMBER_NAME_TO_JSON_VAL_,JSONCONS_NARGS Seq) Seq)
#define NEM_JSON_N_MEMBER_NAME_TO_JSON_VAL_2(Member, Name) \
  {ajson.try_emplace(Name, aval.Member);} \
else \
  {json_traits_helper<Json>::set_optional_json_member(Name, aval.Member, ajson);}
#define NEM_JSON_N_MEMBER_NAME_TO_JSON_VAL_3(Member, Name, Mode) NEM_JSON_N_MEMBER_NAME_TO_JSON_VAL_2(Member, Name)
#define NEM_JSON_N_MEMBER_NAME_TO_JSON_VAL_4(Member, Name, Mode, Match) NEM_JSON_N_MEMBER_NAME_TO_JSON_VAL_6(Member, Name, Mode, Match,,)
#define NEM_JSON_N_MEMBER_NAME_TO_JSON_VAL_5(Member, Name, Mode, Match, Into) NEM_JSON_N_MEMBER_NAME_TO_JSON_VAL_6(Member, Name, Mode, Match, Into, )
#define NEM_JSON_N_MEMBER_NAME_TO_JSON_VAL_6(Member, Name, Mode, Match, Into, From) \
  {Mode ## _TO_JSON_MODE(ajson.try_emplace(Name, Into(aval.Member));)} \
else \
  {Mode ## _TO_JSON_MODE(json_traits_helper<Json>::set_optional_json_member(Name, Into(aval.Member), ajson);)}

#define NEM_JSON_N_MEMBER_NAME_TRAITS_VAL(ValueType,NumMandatoryParams, ...)  \
    JSONCONS_MEMBER_NAME_TRAITS_BASE(JSONCONS_N_MEMBER_NAME_AS, NEM_JSON_N_MEMBER_NAME_TO_JSON_VAL, 0, ValueType,NumMandatoryParams,NumMandatoryParams, __VA_ARGS__) \
    namespace jsoncons { template <> struct is_json_type_traits_declared<ValueType> : public std::true_type {}; } \
  /**/

#define NEM_JSON_CHECK_KEY_ONLY(x) true

#endif  // NEMOSYS_NEMJSONMACROS_H_
