#ifndef LANDAU_LIFSHITZ_HPP
#define LANDAU_LIFSHITZ_HPP

#include <aiwlib/zcube>
#include <aiwlib/gauss>
#include <aiwlib/sphere>
using namespace aiw;

// inline double lang2(double p){ return fabs(p>1e-6)? 1/tanh(p)-1/p: p/3; }
//------------------------------------------------------------------------------
struct Link{
	Ind<3> off; int cell;
	Link(int dx, int dy, int dz, int cell_){ off[0] = dx+1; off[1] = dy+1; off[2] = dz+1; cell = cell_; }
};
//------------------------------------------------------------------------------
extern const char *mode;
// corr это корреляции между крайними атомами в трехчастичной ф.р., при этом расчет ведется относительно атома j

#ifdef CCC
const int cell_sz = 1, nb_sz = 6, corr2_types = 2, corr2_sz[2] = {3,12};
const Link nb_pos[1][6] = {{ Link(-1,0,0, 0), Link(1,0,0, 0), Link(0,-1,0, 0), Link(0,1,0, 0), Link(0,0,-1, 0), Link(0,0,1, 0) }};
const Vecf<3> coord[1] = {vecf(0.f, 0.f, 0.f)};
#   ifdef LL3_CPP
const Ind<2> corr2_0[3] = {ind(0,1), ind(2,3), ind(4,5)};
const Ind<2> corr2_1[12] = {ind(0,2), ind(0,3), ind(0,4), ind(0,5), ind(1,2), ind(1,3), ind(1,4), ind(1,5), ind(2,4), ind(2,5), ind(3,4), ind(3,5)};
const Ind<2> *corr2_table[2] = {corr2_0, corr2_1};
#   endif
#endif
#ifdef VCC
const int cell_sz = 2, nb_sz = 8, corr2_types = 3, corr2_sz[3] = {4,12,12};
const Link nb_pos[2][8] = {{ Link(0,0,0, 1), Link(-1,0,0, 1), Link(0,-1,0, 1), Link(0,0,-1, 1),
							 Link(0,-1,-1, 1), Link(-1,0,-1, 1), Link(-1,-1,0, 1), Link(-1,-1,-1, 1) },
						   { Link(0,0,0, 0), Link(1,0,0, 0), Link(0,1,0, 0), Link(0,0,1, 0),
							 Link(0,1,1, 0), Link(1,0,1, 0), Link(1,1,0, 0), Link(1,1,1, 0) }};
const Vecf<3> coord[2] = {vecf(0.f, 0.f, 0.f), vecf(0.5f, 0.5f, 0.5f)};
#   ifdef LL3_CPP
// autogenerage by make-VCC.py
const Ind<2> corr2_0[4] = {ind(0,7), ind(1,4), ind(2,5), ind(3,6)};
const Ind<2> corr2_1[12] = {ind(0,4), ind(0,5), ind(0,6), ind(1,2), ind(1,3), ind(1,7), ind(2,3), ind(2,7), ind(3,7), ind(4,5), ind(4,6), ind(5,6)};
const Ind<2> corr2_2[12] = {ind(0,1), ind(0,2), ind(0,3), ind(1,5), ind(1,6), ind(2,4), ind(2,6), ind(3,4), ind(3,5), ind(4,7), ind(5,7), ind(6,7)};
const Ind<2> *corr2_table[6] = {corr2_0, corr2_0, corr2_1, corr2_1, corr2_2, corr2_2};
#   endif
#endif
#ifdef FCC
const int cell_sz = 4, nb_sz = 12, corr2_types = 3, corr2_sz[3] = {6, 24, 12};
// autogenerage by make-FCC.py
const Link nb_pos[4][12] = { { Link(-1,-1,0, 3), Link(-1,0,-1, 2), Link(-1,0,0, 2), Link(-1,0,0, 3), Link(0,-1,-1, 1), Link(0,-1,0, 1),
							   Link(0,-1,0, 3), Link(0,0,-1, 1), Link(0,0,-1, 2), Link(0,0,0, 1), Link(0,0,0, 2), Link(0,0,0, 3) },
							 { Link(-1,0,0, 2), Link(-1,0,0, 3), Link(-1,0,1, 3), Link(-1,1,0, 2), Link(0,0,0, 0), Link(0,0,0, 2),
							   Link(0,0,0, 3), Link(0,0,1, 0), Link(0,0,1, 3), Link(0,1,0, 0), Link(0,1,0, 2), Link(0,1,1, 0) },
							 { Link(0,-1,0, 1), Link(0,-1,0, 3), Link(0,-1,1, 3), Link(0,0,0, 0), Link(0,0,0, 1), Link(0,0,0, 3),
							   Link(0,0,1, 0), Link(0,0,1, 3), Link(1,-1,0, 1), Link(1,0,0, 0), Link(1,0,0, 1), Link(1,0,1, 0) },
							 { Link(0,0,-1, 1), Link(0,0,-1, 2), Link(0,0,0, 0), Link(0,0,0, 1), Link(0,0,0, 2), Link(0,1,-1, 2),
							   Link(0,1,0, 0), Link(0,1,0, 2), Link(1,0,-1, 1), Link(1,0,0, 0), Link(1,0,0, 1), Link(1,1,0, 0) }
};
const Vecf<3> coord[4] = {vecf(0.f, 0.f, 0.f), vecf(0.5f, 0.5f, 0.f), vecf(0.5f, 0.f, 0.5f), vecf(0.f, 0.5f, 0.5f)};
#   ifdef LL3_CPP
const Ind<2> corr2_00[6] = { ind(0,11), ind(1,10), ind(2,8), ind(3,6), ind(4,9), ind(5,7) }; // 6
const Ind<2> corr2_01[6] = { ind(0,10), ind(1,8), ind(2,6), ind(3,5), ind(4,11), ind(7,9) }; // 6
const Ind<2> corr2_02[6] = { ind(0,10), ind(1,7), ind(2,5), ind(3,11), ind(4,8), ind(6,9) }; // 6
const Ind<2> corr2_03[6] = { ind(0,10), ind(1,7), ind(2,11), ind(3,8), ind(4,5), ind(6,9) }; // 6

const Ind<2> corr2_10[24] = { ind(0,7), ind(0,8), ind(0,9), ind(0,10), ind(1,5), ind(1,6), ind(1,9), ind(1,11), ind(2,4), ind(2,6), ind(2,7), ind(2,11), ind(3,4), ind(3,5), ind(3,8), ind(3,10), ind(4,10), ind(4,11), ind(5,8), ind(5,11), ind(6,7), ind(6,9), ind(7,10), ind(8,9) }; // 24
const Ind<2> corr2_11[24] = { ind(0,6), ind(0,8), ind(0,9), ind(0,11), ind(1,5), ind(1,7), ind(1,10), ind(1,11), ind(2,4), ind(2,5), ind(2,9), ind(2,10), ind(3,4), ind(3,6), ind(3,7), ind(3,8), ind(4,8), ind(4,10), ind(5,9), ind(5,11), ind(6,7), ind(6,11), ind(7,10), ind(8,9) }; // 24
const Ind<2> corr2_12[24] = { ind(0,5), ind(0,7), ind(0,9), ind(0,11), ind(1,4), ind(1,6), ind(1,10), ind(1,11), ind(2,3), ind(2,4), ind(2,9), ind(2,10), ind(3,7), ind(3,8), ind(3,10), ind(4,9), ind(4,11), ind(5,6), ind(5,8), ind(5,11), ind(6,8), ind(6,10), ind(7,8), ind(7,9) }; // 24
const Ind<2> corr2_13[24] = { ind(0,4), ind(0,7), ind(0,9), ind(0,11), ind(1,3), ind(1,6), ind(1,10), ind(1,11), ind(2,5), ind(2,7), ind(2,8), ind(2,10), ind(3,5), ind(3,9), ind(3,11), ind(4,6), ind(4,8), ind(4,11), ind(5,9), ind(5,10), ind(6,8), ind(6,10), ind(7,8), ind(7,9) }; // 24

const Ind<2> corr2_20[12] = { ind(0,3), ind(0,6), ind(1,2), ind(1,8), ind(2,10), ind(3,11), ind(4,5), ind(4,7), ind(5,9), ind(6,11), ind(7,9), ind(8,10) }; // 12
const Ind<2> corr2_21[12] = { ind(0,3), ind(0,5), ind(1,2), ind(1,6), ind(2,8), ind(3,10), ind(4,7), ind(4,9), ind(5,10), ind(6,8), ind(7,11), ind(9,11) }; // 12
const Ind<2> corr2_22[12] = { ind(0,4), ind(0,8), ind(1,2), ind(1,5), ind(2,7), ind(3,6), ind(3,9), ind(4,10), ind(5,7), ind(6,11), ind(8,10), ind(9,11) }; // 12
const Ind<2> corr2_23[12] = { ind(0,3), ind(0,8), ind(1,4), ind(1,5), ind(2,6), ind(2,9), ind(3,10), ind(4,7), ind(5,7), ind(6,11), ind(8,10), ind(9,11) }; // 12

const Ind<2> *corr2_table[12] = { corr2_00, corr2_01, corr2_02, corr2_03,
								  corr2_10, corr2_11, corr2_12, corr2_13,
								  corr2_20, corr2_21, corr2_22, corr2_23 };
#   endif
#endif
//------------------------------------------------------------------------------
struct Cell{
	Vecf<3> m[cell_sz];
	void init(const Vecf<3> &m0){ for(int i=0; i<cell_sz; i++) m[i] = m0; }
};
//------------------------------------------------------------------------------
class Model{
	File ftm, ftvals; // эволюция одного магнитного момента и tvals
	ZCube<Cell, 3> data[4];

	inline Vecf<3> Hexch(int cube, size_t i, const ZCubeNb<3> &nb, int k) const {
		Vecf<3> H;
		for(const auto &l: nb_pos[k]){ H += data[cube][i+nb[l.off]].m[l.cell];
			WOUT(H, l.off, l.cell, i, i+nb[l.off]);
		}
		return H*J;
	}	
	void calc_av();  // считаем средние значения W, M, M2
	
	std::string path;
	std::vector<Vec<4> > chain_lambda;
	void calc_chain_lambda();
public:
	bool stoch0 = false;
	bool helic0 = false;
	float helic_z = .1;
	float J, gamma, alpha, dt; double T, K, t; // обменный интеграл, температура, прецессия, диссипация, анизотропия, время и шаг по времени
	// int n_q;   // число энергетических уровней при квантовании m
	// float mu = -8;  // химический потенциал в распределении Бозе
	int mg_split = 1;  // минимально возможная длина магнона по каждой из координат (как L/mg_split)
	float zeta;
	aiw::Vecf<3> Hext, nK, M0 = vecf(0.f, 0.f, 1.f);                 // внешнее поле, направление анизотропии и направление начлаьной намагниченности

	bool calc_eq=false; // включает расчет средних
	int eq_count=0, eq_f_count=0;     // число шагов по расчету средних
	aiw::Vec<4> W, Weq; // W, Wexch, Wext, Wanis
	aiw::Vec<3> M, M2, Meq, M2eq;	
	double Mabs_eq = 0.;
	aiw::Vec<3> corr2, corr2eq, Q, Qeq;
	aiw::Vec<4> eta, eta_eq;  // <eta^n>
	void clean_av_eq();
	
	bool f_use = true;
	aiw::Sphere<float> f, f_eq; 

	int data_rank, f_rank;
	void init(const char *path);

	bool parallel = true;
	void calc(int steps);

	void calc_f(); // считаем ф.р.
	void finish();

	bool calc_cl = false; // calc_chain_lambda
	//	void dump(aiw::IOstream &S);
	//	void load(aiw::IOstream &S);
};
//------------------------------------------------------------------------------
#endif //LANDAU_LIFSHITZ_HPP
