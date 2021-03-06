
// Code generated by stanc v2.28.0
#include <stan/model/model_header.hpp>
namespace nonlinear_model_4_model_namespace {

using stan::io::dump;
using stan::model::assign;
using stan::model::index_uni;
using stan::model::index_max;
using stan::model::index_min;
using stan::model::index_min_max;
using stan::model::index_multi;
using stan::model::index_omni;
using stan::model::model_base_crtp;
using stan::model::rvalue;
using namespace stan::math;


stan::math::profile_map profiles__;
static constexpr std::array<const char*, 41> locations_array__ = 
{" (found before start of program)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 10, column 4 to column 13)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 11, column 4 to column 13)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 12, column 4 to column 21)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 13, column 4 to column 21)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 14, column 4 to column 22)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 15, column 4 to column 22)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 16, column 4 to column 23)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 19, column 4 to column 17)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 20, column 4 to column 17)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 21, column 4 to column 22)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 22, column 4 to column 19)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 23, column 4 to column 19)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 26, column 12 to column 48)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 25, column 8 to line 26, column 48)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 24, column 4 to line 26, column 48)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 34, column 4 to column 26)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 37, column 12 to column 75)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 36, column 8 to line 37, column 75)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 35, column 4 to line 37, column 75)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 29, column 4 to column 29)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 30, column 4 to column 29)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 31, column 4 to column 46)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 2, column 4 to column 11)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 3, column 4 to column 11)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 4, column 11 to column 13)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 4, column 4 to column 17)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 5, column 11 to column 13)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 5, column 15 to column 17)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 5, column 4 to column 21)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 6, column 4 to column 14)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 7, column 18 to column 23)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 7, column 4 to column 25)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 12, column 11 to column 13)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 13, column 11 to column 13)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 19, column 11 to column 13)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 20, column 11 to column 13)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 21, column 11 to column 13)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 21, column 15 to column 17)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 34, column 15 to column 17)",
 " (in '/Users/tomoyauchiyama/code/STAN/nonlinear_model_4.stan', line 34, column 19 to column 24)"};



class nonlinear_model_4_model final : public model_base_crtp<nonlinear_model_4_model> {

 private:
  int Np;
  int Nt;
  Eigen::Matrix<double, -1, 1> T__;
  Eigen::Matrix<double, -1, -1> Y__;
  int new_T;
  std::vector<double> new_Time; 
  Eigen::Map<Eigen::Matrix<double, -1, 1>> T{nullptr, 0};
  Eigen::Map<Eigen::Matrix<double, -1, -1>> Y{nullptr, 0, 0};
 
 public:
  ~nonlinear_model_4_model() { }
  
  inline std::string model_name() const final { return "nonlinear_model_4_model"; }

  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.28.0", "stancflags = "};
  }
  
  
  nonlinear_model_4_model(stan::io::var_context& context__,
                          unsigned int random_seed__ = 0,
                          std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double ;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static constexpr const char* function__ = "nonlinear_model_4_model_namespace::nonlinear_model_4_model";
    (void) function__;  // suppress unused var warning
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    try {
      int pos__;
      pos__ = 1;
      current_statement__ = 23;
      context__.validate_dims("data initialization","Np","int",
           std::vector<size_t>{});
      Np = std::numeric_limits<int>::min();
      
      current_statement__ = 23;
      Np = context__.vals_i("Np")[(1 - 1)];
      current_statement__ = 24;
      context__.validate_dims("data initialization","Nt","int",
           std::vector<size_t>{});
      Nt = std::numeric_limits<int>::min();
      
      current_statement__ = 24;
      Nt = context__.vals_i("Nt")[(1 - 1)];
      current_statement__ = 25;
      validate_non_negative_index("T", "Nt", Nt);
      current_statement__ = 26;
      context__.validate_dims("data initialization","T","double",
           std::vector<size_t>{static_cast<size_t>(Nt)});
      T__ = Eigen::Matrix<double, -1, 1>(Nt);
      new (&T) Eigen::Map<Eigen::Matrix<double, -1, 1>>(T__.data(), Nt);
      
      {
        std::vector<local_scalar_t__> T_flat__;
        current_statement__ = 26;
        T_flat__ = context__.vals_r("T");
        current_statement__ = 26;
        pos__ = 1;
        current_statement__ = 26;
        for (int sym1__ = 1; sym1__ <= Nt; ++sym1__) {
          current_statement__ = 26;
          assign(T, T_flat__[(pos__ - 1)],
            "assigning variable T", index_uni(sym1__));
          current_statement__ = 26;
          pos__ = (pos__ + 1);
        }
      }
      current_statement__ = 27;
      validate_non_negative_index("Y", "Np", Np);
      current_statement__ = 28;
      validate_non_negative_index("Y", "Nt", Nt);
      current_statement__ = 29;
      context__.validate_dims("data initialization","Y","double",
           std::vector<size_t>{static_cast<size_t>(Np),
            static_cast<size_t>(Nt)});
      Y__ = Eigen::Matrix<double, -1, -1>(Np, Nt);
      new (&Y) Eigen::Map<Eigen::Matrix<double, -1, -1>>(Y__.data(), Np, Nt);
      
      {
        std::vector<local_scalar_t__> Y_flat__;
        current_statement__ = 29;
        Y_flat__ = context__.vals_r("Y");
        current_statement__ = 29;
        pos__ = 1;
        current_statement__ = 29;
        for (int sym1__ = 1; sym1__ <= Nt; ++sym1__) {
          current_statement__ = 29;
          for (int sym2__ = 1; sym2__ <= Np; ++sym2__) {
            current_statement__ = 29;
            assign(Y, Y_flat__[(pos__ - 1)],
              "assigning variable Y", index_uni(sym2__), index_uni(sym1__));
            current_statement__ = 29;
            pos__ = (pos__ + 1);
          }
        }
      }
      current_statement__ = 30;
      context__.validate_dims("data initialization","new_T","int",
           std::vector<size_t>{});
      new_T = std::numeric_limits<int>::min();
      
      current_statement__ = 30;
      new_T = context__.vals_i("new_T")[(1 - 1)];
      current_statement__ = 31;
      validate_non_negative_index("new_Time", "new_T", new_T);
      current_statement__ = 32;
      context__.validate_dims("data initialization","new_Time","double",
           std::vector<size_t>{static_cast<size_t>(new_T)});
      new_Time = std::vector<double>(new_T, std::numeric_limits<double>::quiet_NaN());
      
      
      current_statement__ = 32;
      new_Time = context__.vals_r("new_Time");
      current_statement__ = 33;
      validate_non_negative_index("log_a", "Np", Np);
      current_statement__ = 34;
      validate_non_negative_index("log_b", "Np", Np);
      current_statement__ = 35;
      validate_non_negative_index("a", "Np", Np);
      current_statement__ = 36;
      validate_non_negative_index("b", "Np", Np);
      current_statement__ = 37;
      validate_non_negative_index("mu", "Np", Np);
      current_statement__ = 38;
      validate_non_negative_index("mu", "Nt", Nt);
      current_statement__ = 39;
      validate_non_negative_index("new_Y", "Np", Np);
      current_statement__ = 40;
      validate_non_negative_index("new_Y", "new_T", new_T);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    num_params_r__ = 1 + 1 + Np + Np + 1 + 1 + 1;
    
  }
  
  template <bool propto__, bool jacobian__ , typename VecR, typename VecI, 
  stan::require_vector_like_t<VecR>* = nullptr, 
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr> 
  inline stan::scalar_type_t<VecR> log_prob_impl(VecR& params_r__,
                                                 VecI& params_i__,
                                                 std::ostream* pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    static constexpr const char* function__ = "nonlinear_model_4_model_namespace::log_prob";
    (void) function__;  // suppress unused var warning
    
    try {
      local_scalar_t__ m_a;
      current_statement__ = 1;
      m_a = in__.template read<local_scalar_t__>();
      local_scalar_t__ m_b;
      current_statement__ = 2;
      m_b = in__.template read<local_scalar_t__>();
      Eigen::Matrix<local_scalar_t__, -1, 1> log_a;
      current_statement__ = 3;
      log_a = in__.template read<Eigen::Matrix<local_scalar_t__, -1, 1>>(Np);
      Eigen::Matrix<local_scalar_t__, -1, 1> log_b;
      current_statement__ = 4;
      log_b = in__.template read<Eigen::Matrix<local_scalar_t__, -1, 1>>(Np);
      local_scalar_t__ s_a;
      current_statement__ = 5;
      s_a = in__.template read_constrain_lb<local_scalar_t__, jacobian__>(0,
              lp__);
      local_scalar_t__ s_b;
      current_statement__ = 6;
      s_b = in__.template read_constrain_lb<local_scalar_t__, jacobian__>(0,
              lp__);
      local_scalar_t__ s_Y;
      current_statement__ = 7;
      s_Y = in__.template read_constrain_lb<local_scalar_t__, jacobian__>(0,
              lp__);
      Eigen::Matrix<local_scalar_t__, -1, 1> a;
      a = Eigen::Matrix<local_scalar_t__, -1, 1>(Np);
      stan::math::fill(a, DUMMY_VAR__);
      
      Eigen::Matrix<local_scalar_t__, -1, 1> b;
      b = Eigen::Matrix<local_scalar_t__, -1, 1>(Np);
      stan::math::fill(b, DUMMY_VAR__);
      
      Eigen::Matrix<local_scalar_t__, -1, -1> mu;
      mu = Eigen::Matrix<local_scalar_t__, -1, -1>(Np, Nt);
      stan::math::fill(mu, DUMMY_VAR__);
      
      current_statement__ = 11;
      assign(a, stan::math::exp(log_a), "assigning variable a");
      current_statement__ = 12;
      assign(b, stan::math::exp(log_b), "assigning variable b");
      current_statement__ = 15;
      for (int j = 1; j <= Nt; ++j) {
        current_statement__ = 14;
        for (int i = 1; i <= Np; ++i) {
          current_statement__ = 13;
          assign(mu,
            (rvalue(a, "a", index_uni(i)) *
              (1 -
                stan::math::exp(
                  (-rvalue(b, "b", index_uni(i)) *
                    rvalue(T, "T", index_uni(j)))))),
            "assigning variable mu", index_uni(i), index_uni(j));
        }
      }
      {
        current_statement__ = 20;
        lp_accum__.add(normal_lpdf<propto__>(log_a, m_a, s_a));
        current_statement__ = 21;
        lp_accum__.add(normal_lpdf<propto__>(log_b, m_b, s_b));
        current_statement__ = 22;
        lp_accum__.add(
          normal_lpdf<propto__>(to_vector(Y), to_vector(mu), s_Y));
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
    } // log_prob_impl() 
    
  template <typename RNG, typename VecR, typename VecI, typename VecVar, 
  stan::require_vector_like_vt<std::is_floating_point, VecR>* = nullptr, 
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr, 
  stan::require_std_vector_vt<std::is_floating_point, VecVar>* = nullptr> 
  inline void write_array_impl(RNG& base_rng__, VecR& params_r__,
                               VecI& params_i__, VecVar& vars__,
                               const bool emit_transformed_parameters__ = true,
                               const bool emit_generated_quantities__ = true,
                               std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    static constexpr bool propto__ = true;
    (void) propto__;
    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    int current_statement__ = 0; 
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    constexpr bool jacobian__ = false;
    (void) DUMMY_VAR__;  // suppress unused var warning
    static constexpr const char* function__ = "nonlinear_model_4_model_namespace::write_array";
    (void) function__;  // suppress unused var warning
    
    try {
      double m_a;
      current_statement__ = 1;
      m_a = in__.template read<local_scalar_t__>();
      double m_b;
      current_statement__ = 2;
      m_b = in__.template read<local_scalar_t__>();
      Eigen::Matrix<double, -1, 1> log_a;
      current_statement__ = 3;
      log_a = in__.template read<Eigen::Matrix<local_scalar_t__, -1, 1>>(Np);
      Eigen::Matrix<double, -1, 1> log_b;
      current_statement__ = 4;
      log_b = in__.template read<Eigen::Matrix<local_scalar_t__, -1, 1>>(Np);
      double s_a;
      current_statement__ = 5;
      s_a = in__.template read_constrain_lb<local_scalar_t__, jacobian__>(0,
              lp__);
      double s_b;
      current_statement__ = 6;
      s_b = in__.template read_constrain_lb<local_scalar_t__, jacobian__>(0,
              lp__);
      double s_Y;
      current_statement__ = 7;
      s_Y = in__.template read_constrain_lb<local_scalar_t__, jacobian__>(0,
              lp__);
      Eigen::Matrix<double, -1, 1> a;
      a = Eigen::Matrix<double, -1, 1>(Np);
      stan::math::fill(a, std::numeric_limits<double>::quiet_NaN());
      
      Eigen::Matrix<double, -1, 1> b;
      b = Eigen::Matrix<double, -1, 1>(Np);
      stan::math::fill(b, std::numeric_limits<double>::quiet_NaN());
      
      Eigen::Matrix<double, -1, -1> mu;
      mu = Eigen::Matrix<double, -1, -1>(Np, Nt);
      stan::math::fill(mu, std::numeric_limits<double>::quiet_NaN());
      
      out__.write(m_a);
      out__.write(m_b);
      out__.write(log_a);
      out__.write(log_b);
      out__.write(s_a);
      out__.write(s_b);
      out__.write(s_Y);
      if (logical_negation((primitive_value(emit_transformed_parameters__) ||
            primitive_value(emit_generated_quantities__)))) {
        return ;
      } 
      current_statement__ = 11;
      assign(a, stan::math::exp(log_a), "assigning variable a");
      current_statement__ = 12;
      assign(b, stan::math::exp(log_b), "assigning variable b");
      current_statement__ = 15;
      for (int j = 1; j <= Nt; ++j) {
        current_statement__ = 14;
        for (int i = 1; i <= Np; ++i) {
          current_statement__ = 13;
          assign(mu,
            (rvalue(a, "a", index_uni(i)) *
              (1 -
                stan::math::exp(
                  (-rvalue(b, "b", index_uni(i)) *
                    rvalue(T, "T", index_uni(j)))))),
            "assigning variable mu", index_uni(i), index_uni(j));
        }
      }
      if (emit_transformed_parameters__) {
        out__.write(a);
        out__.write(b);
        out__.write(mu);
      } 
      if (logical_negation(emit_generated_quantities__)) {
        return ;
      } 
      std::vector<std::vector<double>> new_Y;
      new_Y = std::vector<std::vector<double>>(Np, std::vector<double>(new_T, std::numeric_limits<double>::quiet_NaN()));
      
      
      current_statement__ = 19;
      for (int i = 1; i <= Np; ++i) {
        current_statement__ = 18;
        for (int j = 1; j <= new_T; ++j) {
          current_statement__ = 17;
          assign(new_Y,
            normal_rng(
              (rvalue(a, "a", index_uni(i)) *
                (1 -
                  stan::math::exp(
                    (-rvalue(b, "b", index_uni(i)) *
                      rvalue(new_Time, "new_Time", index_uni(j)))))), s_Y,
              base_rng__),
            "assigning variable new_Y", index_uni(i), index_uni(j));
        }
      }
      for (int sym1__ = 1; sym1__ <= new_T; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= Np; ++sym2__) {
          out__.write(new_Y[(sym2__ - 1)][(sym1__ - 1)]);
        }
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    } // write_array_impl() 
    
  template <typename VecVar, typename VecI, 
  stan::require_std_vector_t<VecVar>* = nullptr, 
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr> 
  inline void transform_inits_impl(VecVar& params_r__, VecI& params_i__,
                                   VecVar& vars__,
                                   std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    
    try {
      int pos__;
      pos__ = 1;
      local_scalar_t__ m_a;
      m_a = in__.read<local_scalar_t__>();
      out__.write(m_a);
      local_scalar_t__ m_b;
      m_b = in__.read<local_scalar_t__>();
      out__.write(m_b);
      Eigen::Matrix<local_scalar_t__, -1, 1> log_a;
      log_a = Eigen::Matrix<local_scalar_t__, -1, 1>(Np);
      stan::math::fill(log_a, DUMMY_VAR__);
      
      for (int sym1__ = 1; sym1__ <= Np; ++sym1__) {
        assign(log_a, in__.read<local_scalar_t__>(),
          "assigning variable log_a", index_uni(sym1__));
      }
      out__.write(log_a);
      Eigen::Matrix<local_scalar_t__, -1, 1> log_b;
      log_b = Eigen::Matrix<local_scalar_t__, -1, 1>(Np);
      stan::math::fill(log_b, DUMMY_VAR__);
      
      for (int sym1__ = 1; sym1__ <= Np; ++sym1__) {
        assign(log_b, in__.read<local_scalar_t__>(),
          "assigning variable log_b", index_uni(sym1__));
      }
      out__.write(log_b);
      local_scalar_t__ s_a;
      s_a = in__.read<local_scalar_t__>();
      out__.write_free_lb(0, s_a);
      local_scalar_t__ s_b;
      s_b = in__.read<local_scalar_t__>();
      out__.write_free_lb(0, s_b);
      local_scalar_t__ s_Y;
      s_Y = in__.read<local_scalar_t__>();
      out__.write_free_lb(0, s_Y);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    } // transform_inits_impl() 
    
  inline void get_param_names(std::vector<std::string>& names__) const {
    
    names__ = std::vector<std::string>{"m_a", "m_b", "log_a", "log_b", "s_a",
      "s_b", "s_Y", "a", "b", "mu", "new_Y"};
    
    } // get_param_names() 
    
  inline void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{},
      std::vector<size_t>{}, std::vector<size_t>{static_cast<size_t>(Np)},
      std::vector<size_t>{static_cast<size_t>(Np)}, std::vector<size_t>{
      }, std::vector<size_t>{}, std::vector<size_t>{},
      std::vector<size_t>{static_cast<size_t>(Np)},
      std::vector<size_t>{static_cast<size_t>(Np)},
      std::vector<size_t>{static_cast<size_t>(Np), static_cast<size_t>(Nt)},
      std::vector<size_t>{static_cast<size_t>(Np), static_cast<size_t>(new_T)}};
    
    } // get_dims() 
    
  inline void constrained_param_names(
                                      std::vector<std::string>& param_names__,
                                      bool emit_transformed_parameters__ = true,
                                      bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "m_a");
    param_names__.emplace_back(std::string() + "m_b");
    for (int sym1__ = 1; sym1__ <= Np; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "log_a" + '.' + std::to_string(sym1__));
      } 
    }
    for (int sym1__ = 1; sym1__ <= Np; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "log_b" + '.' + std::to_string(sym1__));
      } 
    }
    param_names__.emplace_back(std::string() + "s_a");
    param_names__.emplace_back(std::string() + "s_b");
    param_names__.emplace_back(std::string() + "s_Y");
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= Np; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "a" + '.' + std::to_string(sym1__));
        } 
      }
      for (int sym1__ = 1; sym1__ <= Np; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "b" + '.' + std::to_string(sym1__));
        } 
      }
      for (int sym1__ = 1; sym1__ <= Nt; ++sym1__) {
        {
          for (int sym2__ = 1; sym2__ <= Np; ++sym2__) {
            {
              param_names__.emplace_back(std::string() + "mu" + '.' + std::to_string(sym2__) + '.' + std::to_string(sym1__));
            } 
          }
        } 
      }
    }
    
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= new_T; ++sym1__) {
        {
          for (int sym2__ = 1; sym2__ <= Np; ++sym2__) {
            {
              param_names__.emplace_back(std::string() + "new_Y" + '.' + std::to_string(sym2__) + '.' + std::to_string(sym1__));
            } 
          }
        } 
      }
    }
    
    } // constrained_param_names() 
    
  inline void unconstrained_param_names(
                                        std::vector<std::string>& param_names__,
                                        bool emit_transformed_parameters__ = true,
                                        bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "m_a");
    param_names__.emplace_back(std::string() + "m_b");
    for (int sym1__ = 1; sym1__ <= Np; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "log_a" + '.' + std::to_string(sym1__));
      } 
    }
    for (int sym1__ = 1; sym1__ <= Np; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "log_b" + '.' + std::to_string(sym1__));
      } 
    }
    param_names__.emplace_back(std::string() + "s_a");
    param_names__.emplace_back(std::string() + "s_b");
    param_names__.emplace_back(std::string() + "s_Y");
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= Np; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "a" + '.' + std::to_string(sym1__));
        } 
      }
      for (int sym1__ = 1; sym1__ <= Np; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "b" + '.' + std::to_string(sym1__));
        } 
      }
      for (int sym1__ = 1; sym1__ <= Nt; ++sym1__) {
        {
          for (int sym2__ = 1; sym2__ <= Np; ++sym2__) {
            {
              param_names__.emplace_back(std::string() + "mu" + '.' + std::to_string(sym2__) + '.' + std::to_string(sym1__));
            } 
          }
        } 
      }
    }
    
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= new_T; ++sym1__) {
        {
          for (int sym2__ = 1; sym2__ <= Np; ++sym2__) {
            {
              param_names__.emplace_back(std::string() + "new_Y" + '.' + std::to_string(sym2__) + '.' + std::to_string(sym1__));
            } 
          }
        } 
      }
    }
    
    } // unconstrained_param_names() 
    
  inline std::string get_constrained_sizedtypes() const {
    
    return std::string("[{\"name\":\"m_a\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"m_b\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"log_a\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(Np) + "},\"block\":\"parameters\"},{\"name\":\"log_b\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(Np) + "},\"block\":\"parameters\"},{\"name\":\"s_a\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"s_b\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"s_Y\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"a\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(Np) + "},\"block\":\"transformed_parameters\"},{\"name\":\"b\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(Np) + "},\"block\":\"transformed_parameters\"},{\"name\":\"mu\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(Np) + ",\"cols\":" + std::to_string(Nt) + "},\"block\":\"transformed_parameters\"},{\"name\":\"new_Y\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(Np) + ",\"element_type\":{\"name\":\"array\",\"length\":" + std::to_string(new_T) + ",\"element_type\":{\"name\":\"real\"}}},\"block\":\"generated_quantities\"}]");
    
    } // get_constrained_sizedtypes() 
    
  inline std::string get_unconstrained_sizedtypes() const {
    
    return std::string("[{\"name\":\"m_a\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"m_b\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"log_a\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(Np) + "},\"block\":\"parameters\"},{\"name\":\"log_b\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(Np) + "},\"block\":\"parameters\"},{\"name\":\"s_a\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"s_b\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"s_Y\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"a\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(Np) + "},\"block\":\"transformed_parameters\"},{\"name\":\"b\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(Np) + "},\"block\":\"transformed_parameters\"},{\"name\":\"mu\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(Np) + ",\"cols\":" + std::to_string(Nt) + "},\"block\":\"transformed_parameters\"},{\"name\":\"new_Y\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(Np) + ",\"element_type\":{\"name\":\"array\",\"length\":" + std::to_string(new_T) + ",\"element_type\":{\"name\":\"real\"}}},\"block\":\"generated_quantities\"}]");
    
    } // get_unconstrained_sizedtypes() 
    
  
    // Begin method overload boilerplate
    template <typename RNG>
    inline void write_array(RNG& base_rng,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                            const bool emit_transformed_parameters = true,
                            const bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      const size_t num_params__ = 
  ((((((1 + 1) + Np) + Np) + 1) + 1) + 1);
      const size_t num_transformed = 
  ((Np + Np) + (Np * Nt));
      const size_t num_gen_quantities = (Np *
                                                                    new_T);
      std::vector<double> vars_vec(num_params__
       + (emit_transformed_parameters * num_transformed)
       + (emit_generated_quantities * num_gen_quantities));
      std::vector<int> params_i;
      write_array_impl(base_rng, params_r, params_i, vars_vec,
          emit_transformed_parameters, emit_generated_quantities, pstream);
      vars = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>>(
        vars_vec.data(), vars_vec.size());
    }

    template <typename RNG>
    inline void write_array(RNG& base_rng, std::vector<double>& params_r,
                            std::vector<int>& params_i,
                            std::vector<double>& vars,
                            bool emit_transformed_parameters = true,
                            bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      const size_t num_params__ = 
  ((((((1 + 1) + Np) + Np) + 1) + 1) + 1);
      const size_t num_transformed = 
  ((Np + Np) + (Np * Nt));
      const size_t num_gen_quantities = (Np *
                                                                    new_T);
      vars.resize(num_params__
        + (emit_transformed_parameters * num_transformed)
        + (emit_generated_quantities * num_gen_quantities));
      write_array_impl(base_rng, params_r, params_i, vars, emit_transformed_parameters, emit_generated_quantities, pstream);
    }

    template <bool propto__, bool jacobian__, typename T_>
    inline T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
                       std::ostream* pstream = nullptr) const {
      Eigen::Matrix<int, -1, 1> params_i;
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }

    template <bool propto__, bool jacobian__, typename T__>
    inline T__ log_prob(std::vector<T__>& params_r,
                        std::vector<int>& params_i,
                        std::ostream* pstream = nullptr) const {
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }


    inline void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream = nullptr) const final {
      std::vector<double> params_r_vec(params_r.size());
      std::vector<int> params_i;
      transform_inits(context, params_i, params_r_vec, pstream);
      params_r = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>>(
        params_r_vec.data(), params_r_vec.size());
    }

  inline void transform_inits(const stan::io::var_context& context,
                              std::vector<int>& params_i,
                              std::vector<double>& vars,
                              std::ostream* pstream__ = nullptr) const {
     constexpr std::array<const char*, 7> names__{"m_a", "m_b", "log_a",
   "log_b", "s_a", "s_b", "s_Y"};
     const std::array<Eigen::Index, 7> num_params__{1, 1, Np, Np, 1, 
   1, 1};
    
     std::vector<double> params_r_flat__(num_params_r__);
     Eigen::Index size_iter__ = 0;
     Eigen::Index flat_iter__ = 0;
     for (auto&& param_name__ : names__) {
       const auto param_vec__ = context.vals_r(param_name__);
       for (Eigen::Index i = 0; i < num_params__[size_iter__]; ++i) {
         params_r_flat__[flat_iter__] = param_vec__[i];
         ++flat_iter__;
       }
       ++size_iter__;
     }
     vars.resize(params_r_flat__.size());
     transform_inits_impl(params_r_flat__, params_i, vars, pstream__);
    } // transform_inits() 
    
};
}

using stan_model = nonlinear_model_4_model_namespace::nonlinear_model_4_model;

#ifndef USING_R

// Boilerplate
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}

stan::math::profile_map& get_stan_profile_data() {
  return nonlinear_model_4_model_namespace::profiles__;
}

#endif


