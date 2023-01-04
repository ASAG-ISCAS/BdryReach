/**
 * @file   polyReachOptions.cpp
 * @brief  polyReachOptions class
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0
 *
 */

/**
 * @author Changyuan Zhao
 * @date  Nov 2021
 * @version 1.1
 *
 */

#include <dynamics/polyReachOptions.h>

namespace reachSolver
{

	template class polyReachOptions<double>;

	template <typename Number>
	polyReachOptions<Number>::polyReachOptions()
	{
		verbose_ = false;
	}

	template <typename Number>
	const double polyReachOptions<Number>::t() const
	{
		return t_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_t(double t)
	{
		t_ = t;
	}

	template <typename Number>
	const double polyReachOptions<Number>::isRV() const
	{
		return isRV_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_isRV(double isRV)
	{
		isRV_ = isRV;
	}

	template <typename Number>
	const double polyReachOptions<Number>::tStart() const
	{
		return tStart_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_tStart(double tStart)
	{
		tStart_ = tStart;
	}

	template <typename Number>
	const double polyReachOptions<Number>::tFinal() const
	{
		return tFinal_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_tFinal(double tFinal)
	{
		tFinal_ = tFinal;
	}

	template <typename Number>
	const double polyReachOptions<Number>::time_step() const
	{
		return time_step_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_time_step(double time_step)
	{
		time_step_ = time_step;
	}

	template <typename Number>
	const polyZonotope<Number> polyReachOptions<Number>::R0() const
	{
		return R0_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_R0(polyZonotope<Number> R0)
	{
		R0_ = R0;
	}

	template <typename Number>
	const Zonotope<Number> polyReachOptions<Number>::U() const
	{
		return U_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_U(Zonotope<Number> U)
	{
		U_ = U;
	}

	template <typename Number>
	const Zonotope<Number> polyReachOptions<Number>::uTrans_lin() const
	{
		return uTrans_lin_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_uTrans_lin(Zonotope<Number> uTrans_lin)
	{
		uTrans_lin_ = uTrans_lin;
	}
	template <typename Number>
	const Vector_t<Number> polyReachOptions<Number>::uTrans_nonlin() const
	{
		return uTrans_nonlin_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_uTrans_nonlin(Vector_t<Number> uTrans_nonlin)
	{
		uTrans_nonlin_ = uTrans_nonlin;
	}

	template <typename Number>
	const Matrix_t<Number> polyReachOptions<Number>::uTrans_Vec() const
	{
		return uTrans_Vec_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_uTrans_Vec(Matrix_t<Number> uTrans_Vec)
	{
		uTrans_Vec_ = uTrans_Vec;
	}

	template <typename Number>
	const Zonotope<Number> polyReachOptions<Number>::Rtrans() const
	{
		return Rtrans_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_Rtrans(Zonotope<Number> Rtrans)
	{
		Rtrans_ = Rtrans;
	}

	template <typename Number>
	const polyZonotope<Number> polyReachOptions<Number>::Rhom() const
	{
		return Rhom_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_Rhom(polyZonotope<Number> Rhom)
	{
		Rhom_ = Rhom;
	}

	template <typename Number>
	const polyZonotope<Number> polyReachOptions<Number>::Rhom_tp() const
	{
		return Rhom_tp_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_Rhom_tp(polyZonotope<Number> Rhom_tp)
	{
		Rhom_tp_ = Rhom_tp;
	}

	template <typename Number>
	const Zonotope<Number> polyReachOptions<Number>::Raux() const
	{
		return Raux_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_Raux(Zonotope<Number> Raux)
	{
		Raux_ = Raux;
	}

	template <typename Number>
	const Zonotope<Number> polyReachOptions<Number>::Rpar() const
	{
		return Rpar_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_Rpar(Zonotope<Number> Rpar)
	{
		Rpar_ = Rpar;
	}

	template <typename Number>
	const size_t polyReachOptions<Number>::taylor_terms() const
	{
		return taylor_terms_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_taylor_terms(size_t taylor_terms)
	{
		taylor_terms_ = taylor_terms;
	}

	template <typename Number>
	const size_t polyReachOptions<Number>::zonotope_order() const
	{
		return zonotope_order_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_zonotope_order(size_t zonotope_order)
	{
		zonotope_order_ = zonotope_order;
	}

	template <typename Number>
	const size_t polyReachOptions<Number>::intermediate_order() const
	{
		return intermediate_order_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_intermediate_order(size_t intermediate_order)
	{
		intermediate_order_ = intermediate_order;
	}

	template <typename Number>
	const size_t polyReachOptions<Number>::error_order() const
	{
		return error_order_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_error_order(size_t error_order)
	{
		error_order_ = error_order;
	}

	template <typename Number>
	const size_t polyReachOptions<Number>::tensor_order() const
	{
		return tensor_order_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_tensor_order(size_t tensor_order)
	{
		tensor_order_ = tensor_order;
	}

	template <typename Number>
	Number* polyReachOptions<Number>::factor() const
	{
		return factor_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_factor(Number* factor)
	{
		factor_ = factor;
	}

	template <typename Number>
	const std::string polyReachOptions<Number>::alg() const
	{
		return alg_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_alg(std::string alg)
	{
		alg_ = alg;
	}

	template <typename Number>
	const bool polyReachOptions<Number>::verbose() const
	{
		return verbose_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_verbose(bool verbose)
	{
		verbose_ = verbose;
	}

	template <typename Number>
	const Vector_t<Number> polyReachOptions<Number>::max_error() const
	{
		return max_error_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_max_error(Vector_t<Number> max_error)
	{
		max_error_ = max_error;
	}

	template <typename Number>
	const int polyReachOptions<Number>::originContained() const
	{
		return originContained_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_originContained(int originContained)
	{
		originContained_ = originContained;
	}

	template <typename Number>
	const double polyReachOptions<Number>::reductionInterval() const
	{
		return reductionInterval_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_reductionInterval(double reductionInterval)
	{
		reductionInterval_ = reductionInterval;
	}

	template <typename Number>
	const int polyReachOptions<Number>::usekrylovError() const
	{
		return usekrylovError_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_usekrylovError(int usekrylovError)
	{
		usekrylovError_ = usekrylovError;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_maxDepGenOrder(int maxDepGenOrder)
	{
		maxDepGenOrder_ = maxDepGenOrder;
	}

	template <typename Number>
	const int polyReachOptions<Number>::maxDepGenOrder() const
	{
		return maxDepGenOrder_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_maxPolyZonoRatio(double maxPolyZonoRatio)
	{
		maxPolyZonoRatio_ = maxPolyZonoRatio;
	}

	template <typename Number>
	const double polyReachOptions<Number>::maxPolyZonoRatio() const
	{
		return maxPolyZonoRatio_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_restructureTechnique(std::string restructureTechnique)
	{
		restructureTechnique_ = restructureTechnique;
	}

	template <typename Number>
	const std::string polyReachOptions<Number>::restructureTechnique() const
	{
		return restructureTechnique_;
	}

	template <typename Number>
	void polyReachOptions<Number>::set_volApproxMethod(std::string volApproxMethod)
	{
		volApproxMethod_ = volApproxMethod;
	}

	template <typename Number>
	const std::string polyReachOptions<Number>::volApproxMethod() const
	{
		return volApproxMethod_;
	}

} // namespace reachSolver