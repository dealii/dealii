// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

#  include <deal.II/base/signaling_nan.h>

#  include <deal.II/differentiation/sd/symengine_number_types.h>
#  include <deal.II/differentiation/sd/symengine_types.h>
#  include <deal.II/differentiation/sd/symengine_utilities.h>

#  include <symengine/complex_double.h>
#  include <symengine/logic.h>
#  include <symengine/number.h>
#  include <symengine/parser.h>
#  include <symengine/real_double.h>
#  include <symengine/symbol.h>
#  include <symengine/symengine_exception.h>

#  include <string>

DEAL_II_NAMESPACE_OPEN

namespace Differentiation
{
  namespace SD
  {
    namespace SE = ::SymEngine;


    /* ---------------------------- Constructors -------------------------- */


    Expression::Expression()
      : expression()
      , hash(0)
      , hash_computed(false)
    {}


    Expression::Expression(const bool value)
      : expression(SE::boolean(value))
      , hash(0)
      , hash_computed(false)
    {}


    Expression::Expression(const SymEngine::integer_class &value)
      : expression(value)
      , hash(0)
      , hash_computed(false)
    {}


    Expression::Expression(const SymEngine::rational_class &value)
      : expression(value)
      , hash(0)
      , hash_computed(false)
    {}


    Expression::Expression(const Expression &condition,
                           const Expression &expression_if_true,
                           const Expression &expression_if_false)
      : hash(0)
      , hash_computed(false)
    {
      Assert(SE::is_a_Boolean(condition.get_value()),
             ExcMessage(
               "The conditional expression must return a boolean type."));

      const SE::RCP<const SE::Boolean> condition_rcp =
        SE::rcp_static_cast<const SE::Boolean>(condition.get_RCP());
      expression =
        SE::piecewise({{expression_if_true.get_RCP(), condition_rcp},
                       {expression_if_false.get_RCP(), SE::boolTrue}});
    }


    Expression::Expression(const std::vector<std::pair<Expression, Expression>>
                             &               condition_expression,
                           const Expression &expression_otherwise)
      : hash(0)
      , hash_computed(false)
    {
      SE::PiecewiseVec piecewise_function;
      piecewise_function.reserve(condition_expression.size() + 1);

      // Add tested conditional entries
      for (const auto &entry : condition_expression)
        {
          Assert(SE::is_a_Boolean(entry.first.get_value()),
                 ExcMessage(
                   "The conditional expression must return a boolean type."));
          piecewise_function.push_back(
            {entry.second.get_RCP(),
             SE::rcp_static_cast<const SE::Boolean>(entry.first.get_RCP())});
        }

      // Add default value
      piecewise_function.push_back(
        {expression_otherwise.get_RCP(), SE::boolTrue});

      // Initialize
      expression = SE::piecewise(std::move(piecewise_function));
    }


    Expression::Expression(const std::vector<std::pair<Expression, Expression>>
                             &condition_expression)
      : hash(0)
      , hash_computed(false)
    {
      // Use the other constructor with a fatal termination point
      // ensuring that an error is thrown if none of the conditions
      // are met.
      *this = Expression(condition_expression,
                         Expression(numbers::signaling_nan<double>()));
    }


    Expression::Expression(const char *symbol)
      : expression(SE::symbol(symbol))
      , hash(0)
      , hash_computed(false)
    {}


    Expression::Expression(const std::string &str,
                           const bool         parse_as_expression)
      : hash(0)
      , hash_computed(false)
    {
      try
        {
          expression = (parse_as_expression ?
                          SE::parse(str) // The string is a symbolic "name"
                          :
                          SE::rcp_static_cast<const SE::Basic>(SE::symbol(
                            str))); // The string is a symbolic "expression"
        }
      catch (...)
        {
          AssertThrow(false, ExcSymEngineParserError(str));
        }
    }


    Expression::Expression(const std::string &         symbol_func,
                           const types::symbol_vector &arguments)
      : expression(SE::function_symbol(
          symbol_func,
          Utilities::convert_expression_vector_to_basic_vector(arguments)))
      , hash(0)
      , hash_computed(false)
    {}


    Expression::Expression(const SymEngine::Expression &rhs)
      : expression(rhs)
      , hash(0)
      , hash_computed(false)
    {}


    Expression::Expression(const SymEngine::RCP<const SymEngine::Basic> &rhs)
      : expression(rhs)
      , hash(0)
      , hash_computed(false)
    {}


    Expression::Expression(SymEngine::RCP<const SymEngine::Basic> &&rhs)
      : expression(rhs)
      , hash(0)
      , hash_computed(false)
    {}


    /* ------------------------------ Utilities ---------------------------- */


    Expression &
    Expression::parse(const std::string &expression)
    {
      *this = Expression(expression, true /*parse_as_expression*/);
      return *this;
    }


    std::ostream &
    Expression::print(std::ostream &os) const
    {
      os << *this;
      return os;
    }


    void
    Expression::compute_hash()
    {
      // This call is expensive, so don't do any work if we don't have to.
      if (is_hashed())
        return;

      std::hash<SymEngine::Basic> hash_fn;
      hash          = hash_fn(get_value());
      hash_computed = true;
    }


    bool
    Expression::is_hashed() const
    {
      // Note: We cannot rely on the a check like this
      //   return !(hash == 0);
      // because some values might legitimately map to zero.
      return hash_computed;
    }


    SymEngine::hash_t
    Expression::get_hash() const
    {
      Assert(is_hashed(), ExcMessage("Object has not yet been hashed."));
      return hash;
    }


    void
    Expression::reset_hash()
    {
      hash          = 0;
      hash_computed = false;
    }


    void
    Expression::save(std::ostream &os) const
    {
      // We write each expression on a new line.
      // Note: SymEngine outputs a non-terminating string
      os << *this;
      os << std::endl;
    }


    void
    Expression::load(std::istream &is)
    {
      // Need to make sure that we read the entire line in,
      // and then subsequently parse it.
      std::string expr;
      std::getline(is, expr);
      Assert(!is.bad(), ExcIO());
      parse(expr);
    }


    /* ------------------------------- Values ----------------------------- */


    const SE::Expression &
    Expression::get_expression() const
    {
      return expression;
    }


    SE::Expression &
    Expression::get_expression()
    {
      return expression;
    }


    const SE::Basic &
    Expression::get_value() const
    {
      return *get_RCP();
    }


    const SE::RCP<const SE::Basic> &
    Expression::get_RCP() const
    {
      return get_expression().get_basic();
    }


    /* --------------------------- Differentiation ------------------------- */


    Expression
    Expression::differentiate(
      const SymEngine::RCP<const SymEngine::Symbol> &symbol) const
    {
      return Expression(SE::diff(get_RCP(), symbol));
    }


    Expression
    Expression::differentiate(
      const SymEngine::RCP<const SymEngine::Basic> &symbol) const
    {
      // Potential symbol
      return Expression(SE::sdiff(get_RCP(), symbol));
    }


    Expression
    Expression::differentiate(const Expression &symbol) const
    {
      return differentiate(symbol.get_RCP());
    }


    /* ------------- Conversion operators ------------------------- */


    Expression::operator const SymEngine::Expression &() const
    {
      return get_expression();
    }


    Expression::operator const SymEngine::RCP<const SymEngine::Basic> &() const
    {
      return get_expression().get_basic();
    }


    /* ------------- Dictionary-based substitution ------------------------- */


    Expression
    Expression::substitute(
      const SymEngine::map_basic_basic &substitution_values) const
    {
      return Expression(get_expression().subs(substitution_values));
    }


    Expression
    Expression::substitute(
      const types::substitution_map &substitution_values) const
    {
      return substitute(
        Utilities::convert_expression_map_to_basic_map(substitution_values));
    }


    Expression
    Expression::substitute(const Expression &symbol,
                           const Expression &value) const
    {
      Assert(SE::is_a<SE::Symbol>(symbol.get_value()),
             ExcMessage(
               "Substitution with a number that does not represent a symbol."));

      types::substitution_map sub_vals;
      sub_vals[symbol] = value;
      return substitute(sub_vals);
    }


    /* -------------------- Math and relational operators ------------------ */


    Expression &
    Expression::operator=(const Expression &rhs)
    {
      if (this != &rhs)
        {
          this->expression = rhs.get_expression();
          if (rhs.is_hashed())
            {
              this->hash          = rhs.get_hash();
              this->hash_computed = true;
            }
          else
            {
              reset_hash();
            }
        }

      return *this;
    }


    Expression &
    Expression::operator=(Expression &&rhs) noexcept
    {
      if (this != &rhs)
        {
          this->expression    = std::move(rhs.expression);
          this->hash          = std::move(rhs.hash);
          this->hash_computed = std::move(rhs.hash_computed);
        }

      return *this;
    }


    Expression
    Expression::operator-() const
    {
      return Expression(-get_expression());
    }


    Expression &
    Expression::operator+=(const Expression &rhs)
    {
      this->expression += rhs.get_expression();
      reset_hash();
      return *this;
    }


    Expression &
    Expression::operator-=(const Expression &rhs)
    {
      this->expression -= rhs.get_expression();
      reset_hash();
      return *this;
    }


    Expression &
    Expression::operator*=(const Expression &rhs)
    {
      this->expression *= rhs.get_expression();
      reset_hash();
      return *this;
    }


    Expression &
    Expression::operator/=(const Expression &rhs)
    {
      this->expression /= rhs.get_expression();
      reset_hash();
      return *this;
    }


    std::ostream &
    operator<<(std::ostream &stream, const Expression &expr)
    {
      stream << expr.get_expression();
      return stream;
    }


    std::istream &
    operator>>(std::istream &stream, Expression &expr)
    {
      std::string str;
      stream >> str;
      expr.parse(str);
      return stream;
    }


    Expression
    operator==(const Expression &lhs, const Expression &rhs)
    {
      return Expression(SE::Eq(lhs.get_RCP(), rhs.get_RCP()));
    }


    Expression
    operator!=(const Expression &lhs, const Expression &rhs)
    {
      return Expression(SE::Ne(lhs.get_RCP(), rhs.get_RCP()));
    }


    Expression
    operator<(const Expression &lhs, const Expression &rhs)
    {
      return Expression(SE::Lt(lhs.get_RCP(), rhs.get_RCP()));
    }


    Expression
    operator>(const Expression &lhs, const Expression &rhs)
    {
      return Expression(SE::Gt(lhs.get_RCP(), rhs.get_RCP()));
    }


    Expression
    operator<=(const Expression &lhs, const Expression &rhs)
    {
      return Expression(SE::Le(lhs.get_RCP(), rhs.get_RCP()));
    }


    Expression
    operator>=(const Expression &lhs, const Expression &rhs)
    {
      return Expression(SE::Ge(lhs.get_RCP(), rhs.get_RCP()));
    }


    Expression
    operator!(const Expression &expression)
    {
      Assert(SE::is_a_Boolean(expression.get_value()),
             ExcMessage("The expression must return a boolean type."));

      const SE::RCP<const SE::Boolean> expression_rcp =
        SE::rcp_static_cast<const SE::Boolean>(expression.get_RCP());

      return Expression(SE::logical_not(expression_rcp));
    }


    Expression
    operator&(const Expression &lhs, const Expression &rhs)
    {
      Assert(SE::is_a_Boolean(lhs.get_value()),
             ExcMessage("The lhs expression must return a boolean type."));
      Assert(SE::is_a_Boolean(rhs.get_value()),
             ExcMessage("The rhs expression must return a boolean type."));

      const SE::RCP<const SE::Boolean> lhs_rcp =
        SE::rcp_static_cast<const SE::Boolean>(lhs.get_RCP());
      const SE::RCP<const SE::Boolean> rhs_rcp =
        SE::rcp_static_cast<const SE::Boolean>(rhs.get_RCP());

      return Expression(SE::logical_and({lhs_rcp, rhs_rcp}));
    }


    Expression
    operator|(const Expression &lhs, const Expression &rhs)
    {
      Assert(SE::is_a_Boolean(lhs.get_value()),
             ExcMessage("The lhs expression must return a boolean type."));
      Assert(SE::is_a_Boolean(rhs.get_value()),
             ExcMessage("The rhs expression must return a boolean type."));

      const SE::RCP<const SE::Boolean> lhs_rcp =
        SE::rcp_static_cast<const SE::Boolean>(lhs.get_RCP());
      const SE::RCP<const SE::Boolean> rhs_rcp =
        SE::rcp_static_cast<const SE::Boolean>(rhs.get_RCP());

      return Expression(SE::logical_or({lhs_rcp, rhs_rcp}));
    }


    Expression
    operator^(const Expression &lhs, const Expression &rhs)
    {
      Assert(SE::is_a_Boolean(lhs.get_value()),
             ExcMessage("The lhs expression must return a boolean type."));
      Assert(SE::is_a_Boolean(rhs.get_value()),
             ExcMessage("The rhs expression must return a boolean type."));

      const SE::RCP<const SE::Boolean> lhs_rcp =
        SE::rcp_static_cast<const SE::Boolean>(lhs.get_RCP());
      const SE::RCP<const SE::Boolean> rhs_rcp =
        SE::rcp_static_cast<const SE::Boolean>(rhs.get_RCP());

      return Expression(SE::logical_xor({lhs_rcp, rhs_rcp}));
    }


    Expression
    operator&&(const Expression &lhs, const Expression &rhs)
    {
      return lhs & rhs;
    }


    Expression
    operator||(const Expression &lhs, const Expression &rhs)
    {
      return lhs | rhs;
    }


    Expression
    operator+(Expression lhs, const Expression &rhs)
    {
      lhs += rhs;
      return lhs;
    }


    Expression
    operator-(Expression lhs, const Expression &rhs)
    {
      lhs -= rhs;
      return lhs;
    }


    Expression
    operator*(Expression lhs, const Expression &rhs)
    {
      lhs *= rhs;
      return lhs;
    }


    Expression
    operator/(Expression lhs, const Expression &rhs)
    {
      lhs /= rhs;
      return lhs;
    }


  } // namespace SD
} // namespace Differentiation


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE
