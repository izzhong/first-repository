//@header	:	gr_precision.h
//@abstract	:	1.@typedef Real as doulbe/float for need.
//				2.Define some real @constant.
//				3.Define some real deal @function.
//@author	:	zhong
//@date		:	2018/8
//@version	:	bate1.0

#ifndef _GRID_PRECISION_H_
#define _GRID_PRECISION_H_

#include<cfloat>

namespace grid 
{
	//@typedef-->

	//mabye some machines can not support double.
	//so you can change this typedef to [float] if need.
#define REAL_DOUBLE
#ifdef REAL_DOUBLE
	using Real = double;
	constexpr double REAL_MAX = DBL_MAX;
#endif
#ifdef REAL_FLOAT
	using Real float;
	constexpr float REAL_MAX = FLT_MAX;
#endif

	//@constant-->Mathematical constant

	//if the fabs of value < the precision
	//the number will be seen as 0.0
	constexpr Real REAL_PRECISION			=	1e-6;

	//Circumference
	constexpr Real MATH_CONSTANT_PI			=	3.14159265358979323;

	//Reciprocal of Circumference
	constexpr Real MATH_CONSTATN_INV_PI		=	0.31830988618379067;

	//Natural Constant
	constexpr Real MATH_CONSTANT_E			=	2.718281828459045235;

	//@constant-->Physical constant

	//Fundamental physical unit
	//[m]	:	meter
	//[kg]	:	kilogram
	//[s]	:	seconds
	//[A]	:	Ampere
	//[K]	:	kelvin
	//[mol] :	mole
	//[cd]	:	candela
	//Derived unit
	//[N]	:	newton
	//[J]	:	joule
	//[C]	:	coulomb
	//[F]	:	farad
	//[H]	:	henry
	//[T]	:	tesla
	
	//Vacuum speed of light(��չ���)
	//unit	:	m/s 
	constexpr Real PHYCICAL_CONSTANT_c		=	299792458.0;

	//Constant of gravitation(ţ����������)
	//unit	:	m^3/(kg*s^2) 
	constexpr Real PHYCICAL_CONSTANT_G		=	6.67428e-11;

	//Avogadro constant(����٤���޳���)
	//unit	:	mol^-1
	constexpr Real PHYCICAL_CONSTANT_NA		=	6.02214179e23;

	//Universal molar gas constant(����Ħ�����峣��)
	//Unit	:	J/(mol*K)
	constexpr Real PHYCICAL_CONSTANT_R		=	8.314472;
	
	//Boltzman constant(������������)
	//Unit	:	J/K
	//Tips	:	K = R / NA
	constexpr Real PHYCICAL_CONSTANT_K		=	1.3806504e-23;

	//Ideal gas molar volume(��������Ħ�����)
	//Unit	:	m^3
	constexpr Real PHYCICAL_CONSTANT_Vm		=	22.413996e-3;

	//Elementary charge(Ԫ�����)
	//Unit	:	C
	constexpr Real PHYCICAL_CONSTANT_e		=	1.602176487e-19;

	//Atomic mass(ԭ������)
	//Unit	:	kg
	constexpr Real PHYCICAL_CONSTANT_mu		=	1.660538782e-27;

	//Electronic mass(��������)
	//Unit	:	kg
	constexpr Real PHYCICAL_CONSTANT_me		=	9.10938215e-31;

	//Electron charge mass ratio(���Ӻ��ʱ�)
	//Unit	:	C/kg
	//Tips	:	= -e/me
	constexpr Real PHYCICAL_CONSTANT_e_me	=	-1.758820150e11;

	//Proton mass(��������)
	//Unit	:	kg
	constexpr Real PHYCICAL_CONSTANT_mp		=	1.672621637e-27;

	//Neutron mass(��������)
	//Unit	:	kg
	constexpr Real PHYCICAL_CONSTANT_mn		=	1.674927211e-27;

	//Faraday constant(�����ڳ���)
	//Unit	:	C/mol
	//Tips:	:	F = NA * e;
	constexpr Real PHYCICAL_CONSTANT_F		=	9.64853399e4;

	//Permittivity of vacuum(��յ�����)
	//Unit	:	F/m
	constexpr Real PHYCICAL_CONSTANT_e0		=	8.854187817e-12;

	//Permeability of vacuum(��մŵ���)
	//Unit	:	H/m
	constexpr Real PHYCICAL_CONSTANT_u0		=	1.2566370614e-6;

	//Electron magnetic moment(���Ӵž�)
	//Unit	:	J/T
	constexpr Real PHYCICAL_CONSTANT_ue		=	-9.28476377e-24;

	//Proton moment(���Ӵž�)
	//Unit	:	J/T
	constexpr Real PHYCICAL_CONSTANT_up		=	1.4101606662e-26;

	//Bohr radius(�����뾶)
	//Unit	:	m
	constexpr Real PHYCICAL_CONSTANT_a0		=	5.2917720859e-11;

	//Bohr magneton(��������)
	//Unit	:	J/T
	constexpr Real PHYCICAL_CONSTANT_uB		=	9.27400915e-24;

	//Nuclear magneton(�˴���)
	//Unit	:	J/T
	constexpr Real PHYCICAL_CONSTANT_uN		=	5.05078324e-27;

	//Planck constant(���ʿ˳���)
	//Unit	:	J*s
	constexpr Real PHYCICAL_CONSTANT_h		=	6.62606896e-34;

	//Fine structure constant(��ϸ�ṹ����)
	//Unit	:	1
	constexpr Real PHYCICAL_CONSTANT_a		=	7.2973525396e-3;

	//Rydberg constant(��²�����)
	//Unit	:	m^-1
	constexpr Real PHYCICAL_CONSTANT_R00	=	1.0973731568527e7;

	//Compton wavelength(���նٲ���)
	//Unit	:	m
	constexpr Real PHYCICAL_CONSTANT_Yc		=	2.4263102175e-12;

	//Proton electron mass ratio(���ӵ���������)
	//Unit	:	1
	//Tips	:	= mp / me
	constexpr Real PHYCICAL_CONSTANT_mp_me	=	1836.15267247;

	//Electrostatic constant(����������)
	//Unit	:	N*m^2/C^2
	constexpr Real PHYCICAL_CONSTANT_k		=	9.0e9;

	//@function-->

	//@abstract : check the fabs of value is < precisoon or not
	//@tips		: only use when compare
	inline bool IsZeroReal(Real r)
	{	//remove the precision error.
		return (r < 0.0 ? -r : r) < REAL_PRECISION;
	} 

	inline bool IsZeroReal(Real r, Real precision)
	{
		return (r < 0.0 ? -r : r) < precision;
	}

	//@abstract : you can use it like this : _real(real value)
	//to trans a real value so you can remove the precision error.
	//@tips		: only use when compare
	//@warning	: not use it to change the value 
	//			  only use it when compare or dispaly
	inline Real _real(Real r)
	{
		return IsZeroReal(r) ? 0.0 : r;
	}

	//@abstract : Transform angle to radian
	inline Real AngleToRadian(Real angle) 
	{	//radian = angle * pi / 180
		return angle * MATH_CONSTANT_PI * 0.005555555555555556;
	} 

	//@abstract : Transform radian to angle
	inline Real RadianToAngle(Real radian)
	{	//angle = 180 * radian / pi;
		return 180.0 * radian * MATH_CONSTATN_INV_PI;
	} 

} //namespace grid

#endif //_GRID_PRECISION_H_