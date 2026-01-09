#ifndef LANDAU_LIFSHITZ_HPP
#define LANDAU_LIFSHITZ_HPP

#include <fstream>
#include <aiwlib/zcube>
// #include <aiwlib/mesh_as_zcube>
#include <aiwlib/gauss>
#include <aiwlib/sphere>
#include <aiwlib/fftw>
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
const int corr_direct_sz = 3;
const Ind<3> corr_direct[3] = {ind(1, 0, 0), ind(0, 1, 0), ind(0, 0, 1)};
const Vecf<3> coord[1] = { Vecf<3>() };
#endif  // SC

#ifdef BCC
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
const int corr_direct_sz = 4;
const Ind<3> corr_direct[4] = {ind(1, 1, 1), ind(-1, 1, 1), ind(1, -1, 1), ind(1, 1, -1), };
const Vecf<3> coord[2] = { Vecf<3>(), Vecf<3>(.5f) };
#endif  // BCC

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
const int corr_direct_sz = 6;
const Ind<3> corr_direct[6] = {ind(0, 1, 1), ind(1, 0, 1), ind(1, 1, 0), ind(0, -1, 1), ind(-1, 0, 1), ind(-1, 1, 0)};
const Vecf<3> coord[4] = { vecf(0.f,0.f,0.f), vecf(0.f,.5f,.5f), vecf(.5f,0.f,.5f), vecf(.5f,.5f,.0f) };
#endif  // FCC

//------------------------------------------------------------------------------
struct Cell{
	Vecf<3> m[cell_sz];
	void init(const Vecf<3> &m0){ for(int i=0; i<cell_sz; i++) m[i] = m0; }
};
//------------------------------------------------------------------------------
class Model{
	RandN01<float> randN01;
	std::ofstream fspectrum;
	aiw::FFTW<3, float> fftw;
	std::vector<double> spectr_eq; int spectr_eq_count = 0; float spectr_f_max = 0;
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
	double eta_old, dot_eta = 0.;  // для расчета Ts

	std::vector<float> fz, fz_eq, f_eta, f_eta_eq, f_eta_buf;
	std::vector<int> corr_direct_offs;  // массив смещений для расчета корреляционных функций, размер corr_direct_sz*(1<<data_rank*3)*corr_max;
	std::vector<aiw::Vec<4> > corr, corr_eq, corr_buf;  // корреляционная функция
	aiw::File corr_fout;                                // файл для вывода корреляционной функции

	void open_tvals();   // открывает заново файлы tvals.dat, corr.dat и tm0.dat,  вызывает расчет средних
	void drop_tvals();
	bool init_conditions = false;  // флаг задания н.у.
public:
	float J = 1.f;                          ///< обменный интеграл
	float gamma = 1.f;                      ///< гиромагнитное отношение
	float alpha = .1f;                      ///< коэффициент диссипации
	float dt = .01f;                        ///< шаг по времени
	double T = 1.;                          ///< температура
	double K = 0.;                          ///< анизотропия
	aiw::Vecf<3> Hext;                      ///< внешнее поле
	aiw::Vecf<3> nK = vecf(0.f, 0.f, 1.f);  ///< направление анизотропии

	bool patchT = false;       ///< режим температурной поправки  (может работать только при CALC_Q)
	double T_sc_weight = 0.1;  ///< вес в линейном фильтре при вычислении температурной поповки
	double T_sc;               ///< текущая схемная температура (вычисляется только при CALC_Q)
	double Tsc_eq;             ///< равновесная схемная температура
	
	float cL;   ///< множитель в f_k, зависит от размера области ???
	
	double t = 0.;           ///< текущее время
	bool calc_eq = false;    ///< включает расчет средних
	aiw::Vec<4> W;           ///< энергия: W (полная), Wexch (обменная), Wext (вклад внешнего поля), Wanis (анизотропии)
	aiw::Vec<4> Weq;         ///< равновесная энергия: W (полная), Wexch (обменная), Wext (вклад внешнего поля), Wanis (анизотропии)
	aiw::Vec<3> M;           ///< средняя намагниченность
	aiw::Vec<3> Meq;         ///< равновесная средняя намагниченность
	aiw::Vec<3> M2;          ///< вторые моменты намагниченности по осям
	aiw::Vec<3> M2eq;	     ///< равновесные вторые моменты намагниченности по осям
	double Mabs_eq = 0;      ///< равновесный модуль намагниченности (рекомендуется использовать вместо него Meq.abs())
	aiw::Vec<4> Q;           ///< Q_k=<m_i*(m_j%(m_j%m_k))> по координационным сферам
	aiw::Vec<4> Qeq;         ///< равновесные Q_k=<m_i*(m_j%(m_j%m_k))> по координационным сферам
	aiw::Vec<4> eta;         ///< степени ближнего параметра порядка <eta>, <eta^2>, <eta^3>, <eta^4>
	aiw::Vec<4> eta_eq;      ///< равновесные степени ближнего параметра порядка <eta>, <eta^2>, <eta^3>, <eta^4>
	aiw::Vec<4> eta_k2;      ///< степени параметра порядка <eta>, <eta^2>, <eta^3>, <eta^4> для второй координационной сферы
	aiw::Vec<4> eta_k2_eq;   ///< равновесные степени параметра порядка <eta>, <eta^2>, <eta^3>, <eta^4> для второй координационной сферы
	aiw::Vec<4> eta_k3;      ///< степени параметра порядка <eta>, <eta^2>, <eta^3>, <eta^4> для третьей координационной сферы
	aiw::Vec<4> eta_k3_eq;   ///< равновесные степени параметра порядка <eta>, <eta^2>, <eta^3>, <eta^4> для третьей координационной сферы
	aiw::Vec<4> eta_k4;      ///< степени параметра порядка <eta>, <eta^2>, <eta^3>, <eta^4> для четвертой координационной сферы
	aiw::Vec<4> eta_k4_eq;   ///< равновесные степени параметра порядка <eta>, <eta^2>, <eta^3>, <eta^4> для четвертой координационной сферы
	aiw::Vec<3> PHI;         ///< вклад анизотропии в прецессию M в первом уравнении CMD
	aiw::Vec<3> PHIeq;       ///< равновесный вклад анизотропии в прецессию M в первом уравнении CMD
	aiw::Vec<3> THETA;       ///< вклад анизотропии в диссипацию M в первом уравнении CMD
	aiw::Vec<3> THETAeq;     ///< равновесный вклад анизотропии в диссипацию M в первом уравнении CMD
	aiw::Vec<3> UpsilonM;    ///< обменное поле Upsilon*<m> в первом уравнении CMD
	aiw::Vec<3> UpsilonMeq;  ///< равновесное обменное поле Upsilon*<m> в первом уравнении CMD
	aiw::Vec<6> XI;          ///< компоненты матрицы XI: Xi_xx, Xi_yy, Xi_zz, Xi_xy, Xi_xz, Xi_yz (вклад внешнего поля в диссипацию) в первом уравнении CMD
	aiw::Vec<6> XIeq;        ///< равновесные компоненты матрицы XI: Xi_xx, Xi_yy, Xi_zz, Xi_xy, Xi_xz, Xi_yz (вклад внешнего поля в диссипацию) в первом уравнении CMD
	double Psi;              ///< Psi=-<m0_i*(m_j%(m_j%nK))*(m_j*nK)>: вклад анизотропии во втором уравнении CMD
	double Psi_eq = 0;       ///< равновесный Psi=-<m0_i*(m_j%(m_j%nK))*(m_j*nK)>: вклад анизотропии во втором уравнении CMD
	aiw::Vec<3> MxMeq;       ///< равновесное <Mav * Mav>
	int n_s = 0;             ///< число состояний спина при квантовании
	size_t Natoms;
	
	// double S;                ///< энтропия (приближенный расчет через S1 и S2)
	// double Seq = 0;          ///< равновесное значение энтропии
	
	int Ms_start = -1;  ///< с какого ранга (размера большой ячейки) сохранять Ms, -1 --- Ms не вычисляется
	void dump_Ms_arr();
	
	int f_rank = -1;          ///< ранг разбиение сферической сетки для расчета одночастичной функции распределения f(m)
	aiw::Sphere<float> f;     ///< одночастичная функция распределения f(m)
	aiw::Sphere<float> f_eq;  ///< равновесная одночастичной функции распределения f(m)

	int fz_sz = 0;  ///< размер одночастичной функции распределения f1(m_z), если 0 функция не строится
	void dump_fz(const char *path, bool eq) const;  ///< сброс одночастичной функции распределения f1(m_z) в .dat файл

	int f_eta_sz = 0;  ///< размер двухчастичной функции распределения f2(eta), если 0 функция не строится
	void dump_f_eta(const char *path, bool eq) const;  ///< сброс двухчастичной функции распределения f2(eta) в .dat файл

	int corr_max = -1;  ///< размер корреляционной функции (максимальное удаление между частицами в ячейках)
	aiw::Vec<4> get_corr_eq(int i) const { return corr_eq.at(i); }  // для финальной сериализации, возвращает степени eta_far
	
	void clean_av_eq();

	bool out_tm0 = false;  ///< вкл/выкл вывод эволюции одного атома в файл на каждом шаге
	
	int data_rank = 5;  ///< определяет размер моделируемой области как 2^data_rank по каждому измерению
	void init(const char *path);

	// все варианты задания н.у. устанавливают t=0 и открывают заново файлы tvals.dat, corr.dat и tm0.dat
	aiw::Vecf<3> M0 = vecf(0.f, 0.f, 1.f);  ///< начальная намагниченность для всех вариантов н.у.
	int helic_n = 1;                        ///< длина спиновой волны при н.у. helic
	void start_gauss();
	void start_helic(); 

	int threads = 0;  ///< задает число тредов, если =0 то берется системное значение
	void calc(int steps);  ///< считает на steps шагов и в конце расситывает средние. При первом вызвое вызывает start_gauss() если не было другого задания н.у.
	void calc_quant(int steps);  ///< считает на steps шагов и в конце расситывает средние. При первом вызвое вызывает start_gauss() если не было другого задания н.у.
	double rt_init = 0;   ///< суммарное время инициализации
	double rt_calc = 0;   ///< суммарное время расчета
	double rt_diagn = 0;  ///< суммарное время расчета диагностики
	
	void finish();  ///< рассчитывает все равновесные значения

	void dump_data(const char *path) const;  ///< сбрасывает состояние (все магнитные моменты) в бинарном формате
    bool load_data(const char *path);        ///< загружает состояние (все магнитные моменты)  в бинарном формате, трактуется как задание н.у.

	void check_rand(int steps, const char *path);  ///< создает файл с приращениями <M>, <eta> от случайного источника по числу шагов, остальные члены отключены.

	void calc_spectrum(const char *path=nullptr);
};
//------------------------------------------------------------------------------
#endif //LANDAU_LIFSHITZ_HPP
