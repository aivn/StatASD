#pragma once

#include <vector>
#include <fstream>

#include <aiwlib/vec>
#include <aiwlib/gauss>
#include <aiwlib/calc_off>

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

	std::ofstream ftvals;
	
	
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
	float n_b;
	
	void init(const char *path);
	void calc(int Nsteps);
	void finish();	

protected:
	aiw::CalcOff<3, int64_t> calc_off;
	
	aiw::Vecf<3> Hexch(int dID, int64_t mID) const {
		aiw::Vecf<3> H;
		if(!dense_mode){
			for(int a=0; a<3; a++) if(links[mID*3+a]) H += J[0]*data[dID][calc_off(mID, a, 1)];
			for(int a=0; a<3; a++){
				int64_t off = calc_off(mID, a, -1);
				if(links[off*3+a]) H += J[0]*data[dID][off];
			}
		} else {
			for(int a=0; a<3; a++) for(int d=-1; d<=1; d+=2) H += J[0]*data[dID][calc_off(mID, a, d)];
			if(dense_mode>=1) for(const aiw::Ind<3> dpos: stencil1) H += J[1]*data[dID][calc_off(mID, dpos)];
			if(dense_mode>=2) for(const aiw::Ind<3> dpos: stencil2) H += J[2]*data[dID][calc_off(mID, dpos)];
		}
		return H;
	}
};
//------------------------------------------------------------------------------
