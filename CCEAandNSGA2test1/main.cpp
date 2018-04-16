#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <time.h>
#include <algorithm>
#define etax 20
#define etam 20
#define EPS 1e-8

using namespace std;

const double LOWBOUND = 0;
const double UPPBOUND = 1;
const int F1_NVAR = 1;
const int F2_NVAR = 29;
const int IN_POP_SIZE = 100; //100
const int OUT_POP_SIZE = 500;
const double P_CROSS = 0.5;
const double P_MUTAT = 0.02;
const int GENE_NUM = 500;

struct F1_INDIVIDUAL
{
	double m_var[F1_NVAR];
	double m_obj;
	int m_np;
	int m_sp;
	int m_fitness;
	void calObj();
};
struct F1_POPULATION
{
	vector<F1_INDIVIDUAL> m_population;

	void initPopulation();
	void geneticOperation();
	void select(int &p1, int &p2);
	void crossover(const F1_INDIVIDUAL &p1, const F1_INDIVIDUAL &p2, F1_INDIVIDUAL &c);
	void mutation(F1_INDIVIDUAL &ind, double rate);
};
struct F2_INDIVIDUAL
{
	double m_var[F2_NVAR];
	double m_obj;
	int m_np;
	int m_sp;
	int m_fitness;
	F1_INDIVIDUAL f1;
	void calObj();
};
struct F2_POPULATION
{
	vector<F2_INDIVIDUAL> m_population;

	void initPopulation();
	void geneticOperation();
	void select(int &p1, int &p2);
	void crossover(const F2_INDIVIDUAL &p1, const F2_INDIVIDUAL &p2, F2_INDIVIDUAL &c);
	void mutation(F2_INDIVIDUAL &ind, double rate);	
};
struct F_INDIVIDUAL
{
	/*F1_INDIVIDUAL f1;
	F2_INDIVIDUAL f2;*/
	//double m_f1_var[F1_NVAR];
	//double m_f2_var[F2_NVAR];
	double m_f1_obj;
	double m_f2_obj;
};
vector<F_INDIVIDUAL> RES_POP;
vector<F_INDIVIDUAL> TMP_POP;
bool F_INDIVIDUAL_cmp(const F_INDIVIDUAL &ind1, const F_INDIVIDUAL &ind2);
double createZeroToOneRand();
bool isCanAddToTmpPop(const F_INDIVIDUAL &ind);
void addToTmpPop(const F_INDIVIDUAL &ind);
bool isCanAddToResPop(const F_INDIVIDUAL &ind);
void addTmpPopToResPop();
void mysort(int *id, double *dis, int size);
int main()
{
	srand(time(NULL));
	F1_POPULATION f1_pop;
	F2_POPULATION f2_pop;
	f1_pop.initPopulation();
	f2_pop.initPopulation();
	
	F_INDIVIDUAL f_ind[IN_POP_SIZE][IN_POP_SIZE];
	F_INDIVIDUAL f1_opt_objs[IN_POP_SIZE];
	F_INDIVIDUAL f2_opt_objs[IN_POP_SIZE];

	while (RES_POP.size() <= OUT_POP_SIZE)
	{
		TMP_POP.clear();
		for (int gene = 0; gene < GENE_NUM; ++gene)
		{
			for (int i = 0; i < IN_POP_SIZE; ++i)
			{
				f1_pop.m_population[i].calObj();
				for (int j = 0; j < IN_POP_SIZE; ++j)
				{
					f2_pop.m_population[j].f1 = f1_pop.m_population[i];
					f2_pop.m_population[j].calObj();
					f_ind[i][j].m_f1_obj = f1_pop.m_population[i].m_obj;
					f_ind[i][j].m_f2_obj = f2_pop.m_population[j].m_obj;
				}
			}
			for (int i = 0; i < IN_POP_SIZE; ++i)
			{
				int tmp = 0;
				for (int j = 0; j < IN_POP_SIZE; ++j)
				{
					if (F_INDIVIDUAL_cmp(f_ind[i][j], f_ind[i][tmp]) == true)
					{
						tmp = j;
					}
				}
				f1_opt_objs[i].m_f1_obj = f_ind[i][tmp].m_f1_obj;
				f1_opt_objs[i].m_f2_obj = f_ind[i][tmp].m_f2_obj;
			}
			for (int j = 0; j < IN_POP_SIZE; ++j)
			{
				int tmp = 0;
				for (int i = 0; i < IN_POP_SIZE; ++i)
				{
					if (F_INDIVIDUAL_cmp(f_ind[i][j], f_ind[tmp][j]) == true)
					{
						tmp = i;
					}
				}
				f2_opt_objs[j].m_f1_obj = f_ind[tmp][j].m_f1_obj;
				f2_opt_objs[j].m_f2_obj = f_ind[tmp][j].m_f2_obj;
			}
			for (int i = 0; i < IN_POP_SIZE; ++i)
			{
				int cnt1 = 0;
				int cnt2 = 0;
				for (int j = 0; j < IN_POP_SIZE; ++j)
				{
					if (F_INDIVIDUAL_cmp(f1_opt_objs[i], f1_opt_objs[j]) == true)
					{
						cnt1++;
						if (isCanAddToTmpPop(f1_opt_objs[i]) == true)
						{
							addToTmpPop(f1_opt_objs[i]);
						}
					}
					if (F_INDIVIDUAL_cmp(f2_opt_objs[i], f2_opt_objs[j]) == true)
					{
						cnt2++;
						if (isCanAddToTmpPop(f2_opt_objs[i]) == true)
						{
							addToTmpPop(f2_opt_objs[i]);
						}
					}
				}
				f1_pop.m_population[i].m_fitness = cnt1;
				f2_pop.m_population[i].m_fitness = cnt2;
			}
			f1_pop.geneticOperation();
			f2_pop.geneticOperation();
		}
		addTmpPopToResPop();
	}
	if (RES_POP.size() > OUT_POP_SIZE)  //选择距离较大的
	{
		double *dis = new double[RES_POP.size()];
		double *obj = new double[RES_POP.size()];
		int *idx = new int[RES_POP.size()];
		int *idd = new int[RES_POP.size()];
		for (int i = 0; i < RES_POP.size(); ++i)
		{
			idx[i] = i;
			dis[i] = 0;
			obj[i] = RES_POP[i].m_f1_obj;
		}
		mysort(idx, obj, RES_POP.size());
		dis[0] = 1.0e+30;
		dis[RES_POP.size() - 1] = 1.0e+30;
		for (int i = 1; i < RES_POP.size() - 1; ++i)
		{
			dis[i] = (RES_POP[idx[i + 1]].m_f1_obj - RES_POP[idx[i - 1]].m_f1_obj) +
				(RES_POP[idx[i + 1]].m_f2_obj - RES_POP[idx[i - 1]].m_f2_obj);
		}
		mysort(idx, dis, RES_POP.size());
		vector<F_INDIVIDUAL> tmpf;
		int index = 0;
		while (tmpf.size() <= OUT_POP_SIZE)
		{
			tmpf.push_back(RES_POP[idx[index++]]);
		}
		RES_POP.clear();
		for (int i = 0; i < OUT_POP_SIZE; ++i)
		{
			RES_POP.push_back(tmpf[i]);
		}
	}
	ofstream out("outdata.txt");
	for (int i = 0; i < RES_POP.size(); ++i)
	{
		out << RES_POP[i].m_f1_obj << "\t" << RES_POP[i].m_f2_obj << endl;
	}
	return 0;
}

union FloatRand
{
	struct
	{
		unsigned long Frac : 23;
		unsigned long Exp : 8;
		unsigned long Signed : 1;
	} BitArea;
	float Value;
	unsigned long Binary; /* for debug only */
};
double createZeroToOneRand()
{
	union FloatRand r;
	r.BitArea.Signed = 0;
	r.BitArea.Exp = 1;
	r.BitArea.Frac = (rand() * rand()) % 0x800000;
	if (r.BitArea.Frac == 0x7FFFFF)
		r.BitArea.Exp = 0x7D;
	else if (r.BitArea.Frac == 0)
		r.BitArea.Exp = 0x7E;
	else
		r.BitArea.Exp = 0x7E;
	return (double)(r.Value - 0.5)*2.0;

	//return rand() % 1000 / 1000.0;
}

void F1_INDIVIDUAL::calObj()
{
	for (int i = 0; i < F1_NVAR; ++i)
	{
		m_obj = m_var[i];
	}
}

void F2_INDIVIDUAL::calObj()
{
	double g = 0;
	for (int i = 0; i < F2_NVAR; ++i)
	{
		g = g + m_var[i];
	}
	g = g / F2_NVAR;
	g = 1 + 9 * g;
	f1.calObj();
	double h = 1 - sqrt(f1.m_obj / g);
	m_obj = g * h;
}

bool F_INDIVIDUAL_cmp(const F_INDIVIDUAL &ind1, const F_INDIVIDUAL &ind2)  //ind1 支配 ind2 ?
{
	if (ind1.m_f1_obj > ind2.m_f1_obj || ind1.m_f2_obj > ind2.m_f2_obj)
	{
		return false;
	}
	if (ind1.m_f1_obj == ind2.m_f1_obj && ind1.m_f2_obj == ind2.m_f2_obj)
	{
		return false;
	}
	return true;
}
bool isCanAddToTmpPop(const F_INDIVIDUAL &ind)
{
	int cnt = 0;
	for (int i = 0; i < TMP_POP.size(); ++i)
	{
		if (F_INDIVIDUAL_cmp(ind, TMP_POP[i]) == true)
		{
			return true;
		}
		if (F_INDIVIDUAL_cmp(TMP_POP[i], ind) == true)
		{
			return false;
		}
		if (TMP_POP[i].m_f1_obj == ind.m_f1_obj && TMP_POP[i].m_f2_obj == ind.m_f2_obj)
		{
			return false;
		}
		cnt++;
	}
	if (cnt == TMP_POP.size()) return true;
}
void addToTmpPop(const F_INDIVIDUAL &ind)
{
	vector<F_INDIVIDUAL> tmp;
	for (int i = 0; i < TMP_POP.size(); ++i)
	{
		if ((F_INDIVIDUAL_cmp(TMP_POP[i], ind) == true) || 
			(F_INDIVIDUAL_cmp(TMP_POP[i], ind) == false && F_INDIVIDUAL_cmp(ind, TMP_POP[i]) == false))
		{
			tmp.push_back(TMP_POP[i]);
		}
	}
	tmp.push_back(ind);
	TMP_POP.clear();
	for (int i = 0; i < tmp.size(); ++i)
	{
		TMP_POP.push_back(tmp[i]);
	}
	tmp.clear();
}

void F1_POPULATION::geneticOperation()
{
	vector<F1_INDIVIDUAL> tmp;
	for (int i = 0; i < IN_POP_SIZE; ++i)
	{
		int index1, index2;
		select(index1, index2);
		F1_INDIVIDUAL child;
		crossover(m_population[index1], m_population[index2], child);
		mutation(child, P_MUTAT);
		tmp.push_back(child);
	}
	m_population.clear();
	for (int i = 0; i < IN_POP_SIZE; ++i)
	{
		m_population.push_back(tmp[i]);
	}
}
void F1_POPULATION::select(int &p1, int &p2)
{
	int t1 = int(createZeroToOneRand() * IN_POP_SIZE);
	int t2 = int(createZeroToOneRand() * IN_POP_SIZE);
	if (m_population[t1].m_fitness > m_population[t2].m_fitness)
		p1 = t1;
	else
		p1 = t2;
	t1 = int(createZeroToOneRand() * IN_POP_SIZE);
	t2 = int(createZeroToOneRand() * IN_POP_SIZE);
	if (m_population[t1].m_fitness > m_population[t2].m_fitness)
		p2 = t1;
	else
		p2 = t2;
}
void F1_POPULATION::crossover(const F1_INDIVIDUAL &p1, const F1_INDIVIDUAL &p2, F1_INDIVIDUAL &c)
{
	double rand;
	double y1, y2, yl, yu;
	double c1, c2;
	double alpha, beta, betaq;
	double eta_c = etax;
	if (createZeroToOneRand() <= 1.0)
	{
		for (int i = 0; i<F1_NVAR; i++)
		{
			if (createZeroToOneRand() <= P_CROSS)
			{
				if (fabs(p1.m_var[i] - p2.m_var[i]) > EPS)
				{
					if (p1.m_var[i] < p2.m_var[i])
					{
						y1 = p1.m_var[i];
						y2 = p2.m_var[i];
					}
					else
					{
						y1 = p2.m_var[i];
						y2 = p1.m_var[i];
					}
					yl = LOWBOUND;
					yu = UPPBOUND;
					rand = createZeroToOneRand();
					beta = 1.0 + (2.0*(y1 - yl) / (y2 - y1));
					alpha = 2.0 - pow(beta, -(eta_c + 1.0));
					if (rand <= (1.0 / alpha))
					{
						betaq = pow((rand*alpha), (1.0 / (eta_c + 1.0)));
					}
					else
					{
						betaq = pow((1.0 / (2.0 - rand*alpha)), (1.0 / (eta_c + 1.0)));
					}
					c1 = 0.5*((y1 + y2) - betaq*(y2 - y1));
					beta = 1.0 + (2.0*(yu - y2) / (y2 - y1));
					alpha = 2.0 - pow(beta, -(eta_c + 1.0));
					if (rand <= (1.0 / alpha))
					{
						betaq = pow((rand*alpha), (1.0 / (eta_c + 1.0)));
					}
					else
					{
						betaq = pow((1.0 / (2.0 - rand*alpha)), (1.0 / (eta_c + 1.0)));
					}
					c2 = 0.5*((y1 + y2) + betaq*(y2 - y1));
					if (c1<yl)
						c1 = yl;
					if (c2<yl)
						c2 = yl;
					if (c1>yu)
						c1 = yu;
					if (c2>yu)
						c2 = yu;
					if (createZeroToOneRand() <= 0.5)
					{
						c.m_var[i] = c2;
					}
					else
					{
						c.m_var[i] = c1;
					}
				}
				else
				{
					c.m_var[i] = p1.m_var[i];
				}
			}
			else
			{
				c.m_var[i] = p1.m_var[i];
			}
		}
	}
	else
	{
		for (int i = 0; i<F1_NVAR; i++)
		{
			c.m_var[i] = p1.m_var[i];
		}
	}
	return;
}
void F1_POPULATION::mutation(F1_INDIVIDUAL &ind, double rate)
{
	double rnd, delta1, delta2, mut_pow, deltaq;
	double y, yl, yu, val, xy;
	double eta_m = etam;
	int id_rnd = int(createZeroToOneRand()*F1_NVAR);
	for (int j = 0; j < F1_NVAR; j++)
	{
		if (createZeroToOneRand() <= rate)
		{
			y = ind.m_var[j];
			yl = LOWBOUND;
			yu = UPPBOUND;
			delta1 = (y - yl) / (yu - yl);
			delta2 = (yu - y) / (yu - yl);
			rnd = createZeroToOneRand();
			mut_pow = 1.0 / (eta_m + 1.0);
			if (rnd <= 0.5)
			{
				xy = 1.0 - delta1;
				val = 2.0*rnd + (1.0 - 2.0*rnd)*(pow(xy, (eta_m + 1.0)));
				deltaq = pow(val, mut_pow) - 1.0;
			}
			else
			{
				xy = 1.0 - delta2;
				val = 2.0*(1.0 - rnd) + 2.0*(rnd - 0.5)*(pow(xy, (eta_m + 1.0)));
				deltaq = 1.0 - (pow(val, mut_pow));
			}
			y = y + deltaq*(yu - yl);
			if (y<yl)
				y = yl;
			if (y>yu)
				y = yu;
			ind.m_var[j] = y;
		}
	}
	return;
}

void F2_POPULATION::geneticOperation()
{
	vector<F2_INDIVIDUAL> tmp;
	for (int i = 0; i < IN_POP_SIZE; ++i)
	{
		int index1, index2;
		select(index1, index2);
		F2_INDIVIDUAL child;
		crossover(m_population[index1], m_population[index2], child);
		mutation(child, P_MUTAT);
		tmp.push_back(child);
	}
	m_population.clear();
	for (int i = 0; i < IN_POP_SIZE; ++i)
	{
		m_population.push_back(tmp[i]);
	}
}
void F2_POPULATION::select(int &p1, int &p2)
{
	int t1 = int(createZeroToOneRand() * IN_POP_SIZE);
	int t2 = int(createZeroToOneRand() * IN_POP_SIZE);
	if (m_population[t1].m_fitness > m_population[t2].m_fitness)
		p1 = t1;
	else
		p1 = t2;
	t1 = int(createZeroToOneRand() * IN_POP_SIZE);
	t2 = int(createZeroToOneRand() * IN_POP_SIZE);
	if (m_population[t1].m_fitness > m_population[t2].m_fitness)
		p2 = t1;
	else
		p2 = t2;
}
void F2_POPULATION::crossover(const F2_INDIVIDUAL &p1, const F2_INDIVIDUAL &p2, F2_INDIVIDUAL &c)
{
	double rand;
	double y1, y2, yl, yu;
	double c1, c2;
	double alpha, beta, betaq;
	double eta_c = etax;
	if (createZeroToOneRand() <= 1.0)
	{
		for (int i = 0; i<F2_NVAR; i++)
		{
			if (createZeroToOneRand() <= P_CROSS)
			{
				if (fabs(p1.m_var[i] - p2.m_var[i]) > EPS)
				{
					if (p1.m_var[i] < p2.m_var[i])
					{
						y1 = p1.m_var[i];
						y2 = p2.m_var[i];
					}
					else
					{
						y1 = p2.m_var[i];
						y2 = p1.m_var[i];
					}
					yl = LOWBOUND;
					yu = UPPBOUND;
					rand = createZeroToOneRand();
					beta = 1.0 + (2.0*(y1 - yl) / (y2 - y1));
					alpha = 2.0 - pow(beta, -(eta_c + 1.0));
					if (rand <= (1.0 / alpha))
					{
						betaq = pow((rand*alpha), (1.0 / (eta_c + 1.0)));
					}
					else
					{
						betaq = pow((1.0 / (2.0 - rand*alpha)), (1.0 / (eta_c + 1.0)));
					}
					c1 = 0.5*((y1 + y2) - betaq*(y2 - y1));
					beta = 1.0 + (2.0*(yu - y2) / (y2 - y1));
					alpha = 2.0 - pow(beta, -(eta_c + 1.0));
					if (rand <= (1.0 / alpha))
					{
						betaq = pow((rand*alpha), (1.0 / (eta_c + 1.0)));
					}
					else
					{
						betaq = pow((1.0 / (2.0 - rand*alpha)), (1.0 / (eta_c + 1.0)));
					}
					c2 = 0.5*((y1 + y2) + betaq*(y2 - y1));
					if (c1<yl)
						c1 = yl;
					if (c2<yl)
						c2 = yl;
					if (c1>yu)
						c1 = yu;
					if (c2>yu)
						c2 = yu;
					if (createZeroToOneRand() <= 0.5)
					{
						c.m_var[i] = c2;
					}
					else
					{
						c.m_var[i] = c1;
					}
				}
				else
				{
					c.m_var[i] = p1.m_var[i];
				}
			}
			else
			{
				c.m_var[i] = p1.m_var[i];
			}
		}
	}
	else
	{
		for (int i = 0; i<F1_NVAR; i++)
		{
			c.m_var[i] = p1.m_var[i];
		}
	}
	return;
}
void F2_POPULATION::mutation(F2_INDIVIDUAL &ind, double rate)
{
	double rnd, delta1, delta2, mut_pow, deltaq;
	double y, yl, yu, val, xy;
	double eta_m = etam;
	int id_rnd = int(createZeroToOneRand()*F1_NVAR);
	for (int j = 0; j < F1_NVAR; j++)
	{
		if (createZeroToOneRand() <= rate)
		{
			y = ind.m_var[j];
			yl = LOWBOUND;
			yu = UPPBOUND;
			delta1 = (y - yl) / (yu - yl);
			delta2 = (yu - y) / (yu - yl);
			rnd = createZeroToOneRand();
			mut_pow = 1.0 / (eta_m + 1.0);
			if (rnd <= 0.5)
			{
				xy = 1.0 - delta1;
				val = 2.0*rnd + (1.0 - 2.0*rnd)*(pow(xy, (eta_m + 1.0)));
				deltaq = pow(val, mut_pow) - 1.0;
			}
			else
			{
				xy = 1.0 - delta2;
				val = 2.0*(1.0 - rnd) + 2.0*(rnd - 0.5)*(pow(xy, (eta_m + 1.0)));
				deltaq = 1.0 - (pow(val, mut_pow));
			}
			y = y + deltaq*(yu - yl);
			if (y<yl)
				y = yl;
			if (y>yu)
				y = yu;
			ind.m_var[j] = y;
		}
	}
	return;
}

bool isCanAddToResPop(const F_INDIVIDUAL &ind)
{
	int cnt = 0;
	for (int i = 0; i < RES_POP.size(); ++i)
	{
		if (F_INDIVIDUAL_cmp(ind, RES_POP[i]) == true)
		{
			return true;
		}
		if (F_INDIVIDUAL_cmp(RES_POP[i], ind) == true)
		{
			return false;
		}
		if (RES_POP[i].m_f1_obj == ind.m_f1_obj && RES_POP[i].m_f2_obj == ind.m_f2_obj)
		{
			return false;
		}
		cnt++;
	}
	if (cnt == RES_POP.size()) return true;
}
void addTmpPopToResPop()
{
	for (int i = 0; i < TMP_POP.size(); ++i)
	{
		if (isCanAddToResPop(TMP_POP[i]) == true)
		{
			vector<F_INDIVIDUAL> tmp;
			for (int n = 0; n < RES_POP.size(); ++n)
			{
				if ((F_INDIVIDUAL_cmp(RES_POP[n], TMP_POP[i]) == true) ||
					(F_INDIVIDUAL_cmp(RES_POP[n], TMP_POP[i]) == false && F_INDIVIDUAL_cmp(TMP_POP[i], RES_POP[n]) == false))
				{
					tmp.push_back(RES_POP[n]);
				}
			}
			tmp.push_back(TMP_POP[i]);
			RES_POP.clear();
			for (int i = 0; i < tmp.size(); ++i)
			{
				RES_POP.push_back(tmp[i]);
			}
		}
	}
}

void F1_POPULATION::initPopulation()
{
	F1_INDIVIDUAL ind;
	for (int i = 0; i < IN_POP_SIZE; ++i)
	{
		for (int j = 0; j < F1_NVAR; ++j)
		{
			ind.m_var[j] = LOWBOUND + createZeroToOneRand()*(UPPBOUND - LOWBOUND);
		}
		//ind.calObj();
		m_population.push_back(ind);
	}
}
void F2_POPULATION::initPopulation()
{
	F2_INDIVIDUAL ind;
	for (int i = 0; i < IN_POP_SIZE; ++i)
	{
		for (int j = 0; j < F2_NVAR; ++j)
		{
			ind.m_var[j] = LOWBOUND + createZeroToOneRand()*(UPPBOUND - LOWBOUND);
		}
		//ind.calObj();
		m_population.push_back(ind);
	}
}

void mysort(int *id, double *dis, int size)   //降序
{
	for (int i = 0; i < size; ++i)
	{
		for (int j = i + 1; j < size; ++j)
		{
			if (dis[i] < dis[j])
			{
				double tmp = dis[i];
				dis[i] = dis[j];
				dis[j] = tmp;
				int tmp2 = id[i];
				id[i] = id[j];
				id[j] = tmp2;
			}
		}
	}
}