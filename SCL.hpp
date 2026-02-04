#pragma once

#include <vector>
#include <fstream>

#include <aiwlib/vec>
#include <aiwlib/gauss>
#include <aiwlib/fftw>
// #include <aiwlib/calc_off>

//------------------------------------------------------------------------------
class Model{
	aiw::RandN01<float> randN01;
	std::vector<aiw::Vecf<3> > data[4];
	std::vector<bool> links;  // для sparse модели, по три связи на атом, содержит информацию о том включена связь или нет
	std::vector<aiw::Ind<3> > stencil1 = {
		aiw::ind(-1,-1, 0), aiw::ind( 1,-1, 0), aiw::ind(-1, 1, 0), aiw::ind( 1, 1, 0),
		aiw::ind(-1, 0,-1), aiw::ind( 1, 0,-1), aiw::ind(-1, 0, 1), aiw::ind( 1, 0, 1),
		aiw::ind( 0,-1,-1), aiw::ind( 0, 1,-1), aiw::ind( 0,-1, 1), aiw::ind( 0, 1, 1)
	};
	std::vector<aiw::Ind<3> > stencil2 = {
		aiw::ind(-1,-1,-1), aiw::ind( 1,-1,-1), aiw::ind(-1, 1,-1), aiw::ind(-1,-1, 1),
		aiw::ind(-1, 1, 1), aiw::ind( 1,-1, 1), aiw::ind( 1, 1,-1), aiw::ind( 1, 1, 1)
	};

	std::ofstream ftvals, fspectrum, fspectr_av;
	aiw::FFTW<3, float> fftw;
	std::vector<double> spectr_eq; int spectr_eq_count = 0; float spectr_f_max = 0;
	std::string path;
	
	int Nth;  // число тредов
	int eq_count = 0;     // число шагов по расчету средних
public:
	aiw::Vecf<3> J = aiw::vecf(1.f, 1.f, 1.f);   ///< обменный интеграл по трем коорд. сферам
	float gamma = 1.f;                      ///< гиромагнитное отношение
	float alpha = .1f;                      ///< коэффициент диссипации
	float dt = .01f;                        ///< шаг по времени
	double T = 1.;                          ///< температура
	double K = 0.;                          ///< анизотропия
	aiw::Vecf<3> Hext;                      ///< внешнее поле
	aiw::Vecf<3> nK = aiw::vecf(0.f, 0.f, 1.f);  ///< направление анизотропии
	double sparse = 1;                      ///< доля отбракакованных связей в sparce модели
	int dense_mode = 0;                     ///< 0 --- sparse, 1 --- доп связи на ребрах куба (2я коорд сфера), 2 --- доп связи в вершинах куба (3я коорд сфера)
	int periodic = 7;
	aiw::Ind<3> mesh_sz;

	double t = 0.;           ///< текущее время
	bool calc_eq = false;    ///< включает расчет средних
	aiw::Vec<4> W;           ///< энергия: W (полная), Wexch (обменная), Wext (вклад внешнего поля), Wanis (анизотропии)
	aiw::Vec<4> Weq;         ///< равновесная энергия: W (полная), Wexch (обменная), Wext (вклад внешнего поля), Wanis (анизотропии)
	aiw::Vec<3> M;           ///< средняя намагниченность
	aiw::Vec<3> Meq;         ///< равновесная средняя намагниченность
	aiw::Vec<3> M2;          ///< вторые моменты намагниченности по осям
	aiw::Vec<3> M2eq;	     ///< равновесные вторые моменты намагниченности по осям
	double Mabs_eq = 0;      ///< равновесный модуль намагниченности (рекомендуется использовать вместо него Meq.abs())
	double eta, eta_eq = 0;
	float n_b;               ///< среднее число соседей у атома
	double lm_eq;            ///< средняя длина магнонов в равновесии
	double nu_eq;            ///< средняя частота магнонов в равновесии
	double nu2_eq;           ///< средний квадрат частоты магнонов в равновесии

	void init(const char *path);
	void calc(int Nsteps);
	void finish();	

	void calc_spectrum(const char *path=nullptr);
	
protected:
	// aiw::CalcOff<3, int64_t> calc_off;
	aiw::Ind<4> mul;
	int64_t calc_off(int64_t off0, int axis, int pm) const {
		int64_t off = off0 + mul[axis]*pm, c0 = off0/mul[axis+1], c = off/mul[axis+1];
		if(off>=0 && c==c0) return off;
		if(periodic&(1<<axis)){
			if(off<0 || c<c0) return off+mul[axis+1];
			return off-mul[axis+1];
		}
		return -1; // out of range
	}
	int64_t calc_off(int64_t off, const aiw::Ind<3> &dpos) const {
		for(int a=0; a<3; a++){
			off = calc_off(off, a, dpos[a]);
			if(off<0) break; 
		}
		return off;
	}

	aiw::Vecf<3> Hexch(int dID, size_t mID){
		aiw::Vecf<3> H;
		if(!dense_mode){
			for(int a=0; a<3; a++) if(links[mID*3+a]){ int64_t off = calc_off(mID, a, 1); if(off>=0) H += J[0]*data[dID][off]; }
			for(int a=0; a<3; a++){
				int64_t off = calc_off(mID, a, -1);
				if(off>=0 && links[off*3+a]) H += J[0]*data[dID][off];
			}
		} else {
			for(int a=0; a<3; a++) for(int d=-1; d<=1; d+=2){ int64_t off = calc_off(mID, a, d); if(off>=0) H += J[0]*data[dID][off]; }
			if(dense_mode>=1) for(const aiw::Ind<3> dpos: stencil1){ int64_t off = calc_off(mID, dpos); if(off>=0) H += J[1]*data[dID][off]; }
			if(dense_mode>=2) for(const aiw::Ind<3> dpos: stencil2){ int64_t off = calc_off(mID, dpos); if(off>=0) H += J[2]*data[dID][off]; }
		}
		return H;
	}
};
//------------------------------------------------------------------------------
