#ifndef LANDAU_LIFSHITZ_HPP
#define LANDAU_LIFSHITZ_HPP

#include <aiwlib/zcube>
#include <aiwlib/gauss>
#include <aiwlib/sphere>
using namespace aiw;  // ???

#ifdef MAGNONS
#include "magnons.hpp"
static const bool magnons = true;
#else  // MAGNONS
static const bool magnons = false;
#endif //MAGNONS

// inline double lang2(double p){ return fabs(p>1e-6)? 1/tanh(p)-1/p: p/3; }
//------------------------------------------------------------------------------
struct Link{
	Ind<3> off; int cell;
	Link(int dx, int dy, int dz, int cell_){ off[0] = dx+1; off[1] = dy+1; off[2] = dz+1; cell = cell_; }
};
//------------------------------------------------------------------------------
extern const char *mode;

#ifdef SC
const int cell_sz = 1, nb_sz = 6, Q_sz = 2;
const Link nb_pos[1][6] = { { Link(0, 0, -1, 0), Link(0, -1, 0, 0), Link(-1, 0, 0, 0), Link(1, 0, 0, 0), Link(0, 1, 0, 0), Link(0, 0, 1, 0) } };
const Ind<2> Qtable[1][15] = {
	{ ind(1, 0), ind(2, 0), ind(2, 1), ind(3, 0), ind(3, 1), ind(4, 0), ind(4, 2), ind(4, 3), ind(5, 1), ind(5, 2), ind(5, 3), ind(5, 4),
	  ind(3, 2), ind(4, 1), ind(5, 0) }
};
const int Qtable_sz[1][2] = { { 12, 3 } };
const Vec<4> Qcoeff = Vec<4>(4., 1., 0., 0.);
#endif  // SC

#ifdef VCC
const int cell_sz = 2, nb_sz = 8, Q_sz = 3;
const Link nb_pos[2][8] = {
	{ Link(-1, -1, -1, 1), Link(0, -1, -1, 1), Link(-1, 0, -1, 1), Link(0, 0, -1, 1), Link(-1, -1, 0, 1), Link(0, -1, 0, 1), Link(-1, 0, 0, 1), Link(0, 0, 0, 1) },
	{ Link( 0,  0,  0, 0), Link(1,  0,  0, 0), Link( 0, 1,  0, 0), Link(1, 1,  0, 0), Link( 0,  0, 1, 0), Link(1,  0, 1, 0), Link( 0, 1, 1, 0), Link(1, 1, 1, 0) }
};
const Ind<2> Qtable[2][28] = {
	{ ind(1, 0), ind(2, 0), ind(3, 1), ind(3, 2), ind(4, 0), ind(5, 1), ind(5, 4), ind(6, 2), ind(6, 4), ind(7, 3), ind(7, 5), ind(7, 6), ind(2, 1), ind(3, 0),
	  ind(4, 1), ind(4, 2), ind(5, 0), ind(5, 3), ind(6, 0), ind(6, 3), ind(6, 5), ind(7, 1), ind(7, 2), ind(7, 4), ind(4, 3), ind(5, 2), ind(6, 1), ind(7, 0) },
	{ ind(1, 0), ind(2, 0), ind(3, 1), ind(3, 2), ind(4, 0), ind(5, 1), ind(5, 4), ind(6, 2), ind(6, 4), ind(7, 3), ind(7, 5), ind(7, 6), ind(2, 1), ind(3, 0),
	  ind(4, 1), ind(4, 2), ind(5, 0), ind(5, 3), ind(6, 0), ind(6, 3), ind(6, 5), ind(7, 1), ind(7, 2), ind(7, 4), ind(4, 3), ind(5, 2), ind(6, 1), ind(7, 0) }
};
const int Qtable_sz[2][3] = { { 12, 12, 4 }, { 12, 12, 4 } };
const Vec<4> Qcoeff = Vec<4>(3., 3., 1., 0.);
#endif  // VCC

#ifdef FCC
const int cell_sz = 4, nb_sz = 12, Q_sz = 4;
const Link nb_pos[4][12] = {
	{ Link(0, -1, -1, 1), Link(-1, 0, -1, 2), Link( 0, 0, -1, 1), Link(0, 0, -1, 2), Link(-1, -1, 0, 3), Link(0, -1, 0, 1),
	  Link(0, -1,  0, 3), Link(-1, 0,  0, 2), Link(-1, 0,  0, 3), Link(0, 0,  0, 1), Link( 0,  0, 0, 2), Link(0,  0, 0, 3) },
	{ Link(-1, 0, 0, 2), Link(-1, 0, 0, 3), Link( 0, 0, 0, 0), Link(0, 0, 0, 2), Link(0, 0, 0, 3), Link(-1, 1, 0, 2),
	  Link( 0, 1, 0, 0), Link( 0, 1, 0, 2), Link(-1, 0, 1, 3), Link(0, 0, 1, 0), Link(0, 0, 1, 3), Link( 0, 1, 1, 0) },
	{ Link(0, -1, 0, 1), Link(0, -1, 0, 3), Link(1, -1, 0, 1), Link(0, 0, 0, 0), Link(0, 0, 0, 1), Link(0, 0, 0, 3),
	  Link(1,  0, 0, 0), Link(1,  0, 0, 1), Link(0, -1, 1, 3), Link(0, 0, 1, 0), Link(0, 0, 1, 3), Link(1, 0, 1, 0) },
	{ Link(0, 0, -1, 1), Link(0, 0, -1, 2), Link(1, 0, -1, 1), Link(0, 1, -1, 2), Link(0, 0, 0, 0), Link(0, 0, 0, 1),
	  Link(0, 0,  0, 2), Link(1, 0,  0, 0), Link(1, 0,  0, 1), Link(0, 1,  0, 0), Link(0, 1, 0, 2), Link(1, 1, 0, 0) }
};
const Ind<2> Qtable[4][66] = {
	{ ind(1, 0), ind(2, 1), ind(3, 0), ind(3, 2), ind(4, 0), ind(4, 1), ind(5, 4), ind(6, 0), ind(6, 3), ind(6, 5), ind(7, 4), ind(7, 5), ind(8, 1), ind(8, 2), ind(8, 7), ind(9, 7), ind(9, 8), ind(10, 5), ind(10, 6), ind(10, 9), ind(11, 2), ind(11, 3), ind(11, 9), ind(11, 10), ind(2, 0), ind(3, 1), ind(5, 0), ind(6, 4), ind(7, 1), ind(8, 4), ind(9, 2), ind(9, 5), ind(10, 3), ind(10, 7), ind(11, 6), ind(11, 8), ind(4, 2), ind(4, 3), ind(5, 1), ind(5, 3), ind(6, 1), ind(6, 2), ind(7, 0), ind(7, 2), ind(7, 6), ind(8, 0), ind(8, 3), ind(8, 5), ind(9, 1), ind(9, 3), ind(9, 4), ind(9, 6), ind(10, 0), ind(10, 2), ind(10, 4), ind(10, 8), ind(11, 0), ind(11, 1), ind(11, 5), ind(11, 7), ind(5, 2), ind(7, 3), ind(8, 6), ind(9, 0), ind(10, 1), ind(11, 4) },
	{ ind(1, 0), ind(2, 0), ind(2, 1), ind(3, 2), ind(4, 2), ind(4, 3), ind(5, 1), ind(6, 1), ind(6, 4), ind(6, 5), ind(7, 4), ind(7, 6), ind(8, 0), ind(8, 5), ind(9, 0), ind(9, 3), ind(9, 8), ind(10, 3), ind(10, 7), ind(10, 9), ind(11, 5), ind(11, 7), ind(11, 8), ind(11, 10), ind(3, 0), ind(4, 1), ind(5, 0), ind(6, 2), ind(7, 3), ind(7, 5), ind(8, 1), ind(9, 2), ind(10, 4), ind(10, 8), ind(11, 6), ind(11, 9), ind(3, 1), ind(4, 0), ind(5, 2), ind(5, 4), ind(6, 0), ind(6, 3), ind(7, 1), ind(7, 2), ind(8, 2), ind(8, 3), ind(8, 6), ind(8, 7), ind(9, 1), ind(9, 4), ind(9, 5), ind(9, 7), ind(10, 0), ind(10, 2), ind(10, 5), ind(10, 6), ind(11, 0), ind(11, 1), ind(11, 3), ind(11, 4), ind(5, 3), ind(7, 0), ind(8, 4), ind(9, 6), ind(10, 1), ind(11, 2) },
	{ ind(1, 0), ind(2, 1), ind(3, 0), ind(3, 1), ind(4, 3), ind(5, 3), ind(5, 4), ind(6, 1), ind(6, 2), ind(6, 5), ind(7, 5), ind(7, 6), ind(8, 0), ind(8, 2), ind(9, 0), ind(9, 4), ind(9, 8), ind(10, 4), ind(10, 7), ind(10, 9), ind(11, 2), ind(11, 7), ind(11, 8), ind(11, 10), ind(2, 0), ind(4, 0), ind(5, 1), ind(6, 3), ind(7, 2), ind(7, 4), ind(8, 1), ind(9, 3), ind(10, 5), ind(10, 8), ind(11, 6), ind(11, 9), ind(3, 2), ind(4, 1), ind(5, 0), ind(5, 2), ind(6, 0), ind(6, 4), ind(7, 1), ind(7, 3), ind(8, 3), ind(8, 4), ind(8, 6), ind(8, 7), ind(9, 1), ind(9, 2), ind(9, 5), ind(9, 7), ind(10, 0), ind(10, 2), ind(10, 3), ind(10, 6), ind(11, 0), ind(11, 1), ind(11, 4), ind(11, 5), ind(4, 2), ind(7, 0), ind(8, 5), ind(9, 6), ind(10, 1), ind(11, 3) },
	{ ind(1, 0), ind(2, 1), ind(3, 0), ind(3, 2), ind(4, 0), ind(4, 1), ind(5, 4), ind(6, 4), ind(6, 5), ind(7, 1), ind(7, 2), ind(7, 6), ind(8, 6), ind(8, 7), ind(9, 0), ind(9, 3), ind(9, 5), ind(10, 5), ind(10, 8), ind(10, 9), ind(11, 2), ind(11, 3), ind(11, 8), ind(11, 10), ind(2, 0), ind(3, 1), ind(5, 0), ind(6, 1), ind(7, 4), ind(8, 2), ind(8, 5), ind(9, 4), ind(10, 3), ind(10, 6), ind(11, 7), ind(11, 9), ind(4, 2), ind(4, 3), ind(5, 1), ind(5, 3), ind(6, 0), ind(6, 2), ind(7, 0), ind(7, 3), ind(7, 5), ind(8, 1), ind(8, 3), ind(8, 4), ind(9, 1), ind(9, 2), ind(9, 6), ind(9, 8), ind(10, 0), ind(10, 2), ind(10, 4), ind(10, 7), ind(11, 0), ind(11, 1), ind(11, 5), ind(11, 6), ind(5, 2), ind(6, 3), ind(8, 0), ind(9, 7), ind(10, 1), ind(11, 4) }
};
const int Qtable_sz[4][4] = { { 24, 12, 24, 6 }, { 24, 12, 24, 6 }, { 24, 12, 24, 6 }, { 24, 12, 24, 6 } };
const Vec<4> Qcoeff = Vec<4>(4., 2., 4., 1.);
#endif  // FCC

//------------------------------------------------------------------------------
struct Cell{
	Vecf<3> m[cell_sz];
	void init(const Vecf<3> &m0){ for(int i=0; i<cell_sz; i++) m[i] = m0; }
};
//------------------------------------------------------------------------------
class Model{
	RandN01<float> randN01;
	aiw::File ftm, ftvals; // эволюция одного магнитного момента и tvals
	ZCube<Cell, 3> data[4];

	inline Vecf<3> Hexch(int cube, size_t i, const ZCubeNb<3> &nb, int k) const {
		Vecf<3> H;
		for(const auto &l: nb_pos[k]){ H += data[cube][i+nb[l.off]].m[l.cell];
			// WOUT(H, l.off, l.cell, i, i+nb[l.off]);
		}
		return H*J;
	}	
	void calc_av();  // считаем средние значения W, M, M2 и тд
	
	std::string path;

#ifdef MAGNONS
	MagnonsStochFieldGen MG;
#endif // MAGNONS

	std::vector<int> f_buf, fz_buf;
	int Nth;  // число тредов
	int eq_count = 0;     // число шагов по расчету средних

	std::vector<aiw::Mesh<aiw::Vecf<3>, 3> > Ms_arr;
	std::vector<aiw::File> fMs;  // файлы для сброса Ms

	std::vector<double> Q_buf, eta_k_buf;
	double eta_old, dot_eta;  // эффективная схемная температура (может работать при CALC_Q)
public:
	bool stoch0 = false;
	bool helic0 = false;
	int helic_n = 1;   // длина спиновой волны
	bool entropy0 = false;  // старт с гауссовым распределением с нулевой энтропией
	float helic_z = .1;
	float J, gamma, alpha, dt; double T, K, t;  // обменный интеграл, температура, прецессия, диссипация, анизотропия, время и шаг по времени

	bool patchT = false;   // режим температурной поправки  (может работать при CALC_Q)
	double T_sc_weight, T_sc, Tsc_eq;
	
	float cL;   // множитель в f_k, зависит от размера области ???
	
	aiw::Vecf<3> Hext, nK, M0 = vecf(0.f, 0.f, 1.f);                 // внешнее поле, направление анизотропии и направление начлаьной намагниченности

	bool calc_eq = false; // включает расчет средних
	aiw::Vec<4> W, Weq; // W, Wexch, Wext, Wanis
	aiw::Vec<3> M, M2, Meq, M2eq;	
	double Mabs_eq = 0.;
	aiw::Vec<4> Q, Qeq;
	aiw::Vec<4> eta, eta_eq, eta_k2, eta_k2_eq, eta_k3, eta_k3_eq, eta_k4, eta_k4_eq;  // <eta^n>
	aiw::Vec<3> PHI, PHIeq, THETA, THETAeq, UpsilonM, UpsilonMeq;
	aiw::Vec<6> XI, XIeq;  
	double Psi, Psi_eq;  

	std::vector<double> Ms, Ms_eq;  // массив средних модулей намагниченности для разных масштабов, [0] - пары ближ. соседей, [1] - удвоенная ячейка, .back() - весь обр.
	void dump_Ms_arr();
	
	aiw::Sphere<float> f, f_eq;

	int fz_sz = 0;  // размер f1(m_z), если 0 не используется
	std::vector<float> fz, fz_eq;
	void dump_fz(const char *path, bool eq) const;

	void clean_av_eq();

	bool out_tm0 = false; // вывод эволюции одного атома
	bool Hinv = false;    // внешнее поле антиколлинеарно <m>
	
	int data_rank, f_rank = -1;
	void init(const char *path);

	bool parallel = true;
	void calc(int steps);

	void finish();
};
//------------------------------------------------------------------------------
#endif //LANDAU_LIFSHITZ_HPP
