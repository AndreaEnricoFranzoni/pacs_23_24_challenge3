#include <muParser.h>
#include <memory>
#include <string>
#include <numbers>

constexpr double pi_ = std::numbers::pi;

class MuparserFun2
{
public:
  MuparserFun2(const MuparserFun2 &m)
    : m_parser(m.m_parser)
  { //the function has two input variables
    m_parser.DefineVar("x", &m_var1);
    m_parser.DefineVar("y", &m_var2);
    m_parser.DefineConst("_pi", pi_);
  };

  MuparserFun2(const std::string &s)
  {
    try
      {
        m_parser.DefineVar("x", &m_var1);
        m_parser.DefineVar("y", &m_var2);
        m_parser.SetExpr(s);
      }
    catch (mu::Parser::exception_type &e)
      {
        std::cerr << e.GetMsg() << std::endl;
      }
  };

  double
  operator()(const double &x, const double &y)
  {
    m_var1 = x;
    m_var2 = y;
    double f_ev = m_parser.Eval();
    return f_ev;
  };

private:
  double     m_var1;
  double     m_var2;
  mu::Parser m_parser;
};
