#include <Rcpp.h>
#include <time.h>
#include <fstream>
using namespace Rcpp;

// Generator of pseudorandom numbers

#define norm 2.328306549295728e-10
#define m1   4294967087.0
#define m2   4294944443.0
#define a12     1403580.0
#define a13n     810728.0
#define a21      527612.0
#define a23n    1370589.0

double  s10, s11, s12, s20, s21, s22;

double MRG32k3a()
{
  long   k;
  double p1, p2;
	/* Component 1 */
	p1 = a12 * s11 - a13n * s10;
	k = p1 / m1;   p1 -= k * m1;   if (p1 < 0.0) p1 += m1;
	s10 = s11;   s11 = s12;   s12 = p1;
	/* Component 2 */
	p2 = a21 * s22 - a23n * s20;
	k =  p2 / m2;  p2 -= k * m2;   if (p2 < 0.0) p2 += m2;
	s20 = s21;   s21 = s22;   s22 = p2;
	/* Combination */
	if (p1 <= p2) return ((p1 - p2 + m1) * norm);
	else return ((p1 - p2) * norm);
}

int get_next()
{
  return floor(64.0 * MRG32k3a());
}

void init_gen()
{
  int p1, p2, p3, p4, p5, p6;
	do
	{
		p1 = rand();
		p2 = rand();
		p3 = rand();
		p4 = rand();
		p5 = rand();
		p6 = rand();
	} while ((p1 == 0) && (p2 == 0) && (p3 == 0) && (p4 == 0) && (p5 == 0) && (p6 == 0));
}

// Declaration of six dimensional qube

class Q6
{
public:

  std::ofstream output_dml; 
  
  Q6() 
  {
    output_dml.open("Founded_DML.txt",std::ios::app); // insert target address
  }
  ~Q6() 
  {
    output_dml.close();
  }
  
  static const int magic_constant = 189;

	bool simuluj();
	inline bool Nastav();
	inline bool Dopocitej_kompl_baze();
	inline bool Dopocitej_pom_baze();
	inline bool Dopocitej_zbyle();

	inline bool mimo_interval(int x) { return ((x < 0) && (x > 63)); }
	inline bool pouzite(int index) { return hodnoty[index]; }

	int vrcholy[64];
	bool hodnoty[64];
	inline void init()
	{
		for (size_t i = 0; i < 64; i++)
		{
			hodnoty[i] = false;
		}
	}
  
  void save_dml()
  {
    for (size_t i = 0; i < 64; i++)
    {
      output_dml << " " << vrcholy[0];
    }
    output_dml << "\n";
  }
};

// Setting of a base of Q_6 (base is A_6), 20 vertices

inline bool Q6::Nastav()
{
  vrcholy[0] = get_next();
	hodnoty[vrcholy[0]] = true;
	vrcholy[1] = get_next();
	if (pouzite(vrcholy[1])) return false;
	hodnoty[vrcholy[1]] = true;
	vrcholy[2] = get_next();
	if (pouzite(vrcholy[2])) return false;
	hodnoty[vrcholy[2]] = true;
	vrcholy[3] = get_next();
	if (pouzite(vrcholy[3])) return false;
	hodnoty[vrcholy[3]] = true;
	vrcholy[4] = get_next();
	if (pouzite(vrcholy[4])) return false;
	hodnoty[vrcholy[4]] = true;
	vrcholy[5] = get_next();
	if (pouzite(vrcholy[5])) return false;
	hodnoty[vrcholy[5]] = true;
	vrcholy[6] = get_next();
	if (pouzite(vrcholy[6])) return false;
	hodnoty[vrcholy[6]] = true;
	vrcholy[7] = get_next();
	if (pouzite(vrcholy[7])) return false;
	hodnoty[vrcholy[7]] = true;
	vrcholy[8] = get_next();
	if (pouzite(vrcholy[8])) return false;
	hodnoty[vrcholy[8]] = true;
	vrcholy[9] = get_next();
	if (pouzite(vrcholy[9])) return false;
	hodnoty[vrcholy[9]] = true;
	vrcholy[10] = get_next();
	if (pouzite(vrcholy[10])) return false;
	hodnoty[vrcholy[10]] = true;
	vrcholy[11] = get_next();
	if (pouzite(vrcholy[11])) return false;
	hodnoty[vrcholy[11]] = true;
	vrcholy[12] = get_next();
	if (pouzite(vrcholy[12])) return false;
	hodnoty[vrcholy[12]] = true;
	vrcholy[13] = get_next();
	if (pouzite(vrcholy[13])) return false;
	hodnoty[vrcholy[13]] = true;
	vrcholy[16] = get_next();
	if (pouzite(vrcholy[16])) return false;
	hodnoty[vrcholy[16]] = true;
	vrcholy[17] = get_next();
	if (pouzite(vrcholy[17])) return false;
	hodnoty[vrcholy[17]] = true;
	vrcholy[18] = get_next();
	if (pouzite(vrcholy[18])) return false;
	hodnoty[vrcholy[18]] = true;
	vrcholy[19] = get_next();
	if (pouzite(vrcholy[19])) return false;
	hodnoty[vrcholy[19]] = true;
	vrcholy[20] = get_next();
	if (pouzite(vrcholy[20])) return false;
	hodnoty[vrcholy[20]] = true;
	vrcholy[21] = get_next();
	if (pouzite(vrcholy[21])) return false;
	hodnoty[vrcholy[21]] = true;
  return true;
}

// Setting of a complement of base, 20 vertices

inline bool Q6::Dopocitej_kompl_baze()
{
  vrcholy[63] = 63 - vrcholy[0];
	if (pouzite(vrcholy[63])) return false;
	hodnoty[vrcholy[63]] = true;
	vrcholy[62] = 63 - vrcholy[1];
	if (pouzite(vrcholy[62])) return false;
	hodnoty[vrcholy[62]] = true;
	vrcholy[61] = 63 - vrcholy[2];
	if (pouzite(vrcholy[61])) return false;
	hodnoty[vrcholy[61]] = true;
	vrcholy[60] = 63 - vrcholy[3];
	if (pouzite(vrcholy[60])) return false;
	hodnoty[vrcholy[60]] = true;
	vrcholy[59] = 63 - vrcholy[4];
	if (pouzite(vrcholy[59])) return false;
	hodnoty[vrcholy[59]] = true;
	vrcholy[58] = 63 - vrcholy[5];
	if (pouzite(vrcholy[58])) return false;
	hodnoty[vrcholy[58]] = true;
	vrcholy[57] = 63 - vrcholy[6];
	if (pouzite(vrcholy[57])) return false;
	hodnoty[vrcholy[57]] = true;
	vrcholy[56] = 63 - vrcholy[7];
	if (pouzite(vrcholy[56])) return false;
	hodnoty[vrcholy[56]] = true;
	vrcholy[55] = 63 - vrcholy[8];
	if (pouzite(vrcholy[55])) return false;
	hodnoty[vrcholy[55]] = true;
	vrcholy[54] = 63 - vrcholy[9];
	if (pouzite(vrcholy[54])) return false;
	hodnoty[vrcholy[54]] = true;
	vrcholy[53] = 63 - vrcholy[10];
	if (pouzite(vrcholy[53])) return false;
	hodnoty[vrcholy[53]] = true;
	vrcholy[52] = 63 - vrcholy[11];
	if (pouzite(vrcholy[52])) return false;
	hodnoty[vrcholy[52]] = true;
	vrcholy[51] = 63 - vrcholy[12];
	if (pouzite(vrcholy[51])) return false;
	hodnoty[vrcholy[51]] = true;
	vrcholy[50] = 63 - vrcholy[13];
	if (pouzite(vrcholy[50])) return false;
	hodnoty[vrcholy[50]] = true;
	vrcholy[47] = 63 - vrcholy[16];
	if (pouzite(vrcholy[47])) return false;
	hodnoty[vrcholy[47]] = true;
	vrcholy[46] = 63 - vrcholy[17];
	if (pouzite(vrcholy[46])) return false;
	hodnoty[vrcholy[46]] = true;
	vrcholy[45] = 63 - vrcholy[18];
	if (pouzite(vrcholy[45])) return false;
	hodnoty[vrcholy[45]] = true;
	vrcholy[44] = 63 - vrcholy[19];
	if (pouzite(vrcholy[44])) return false;
	hodnoty[vrcholy[44]] = true;
	vrcholy[43] = 63 - vrcholy[20];
	if (pouzite(vrcholy[43])) return false;
	hodnoty[vrcholy[43]] = true;
	vrcholy[42] = 63 - vrcholy[21];
	if (pouzite(vrcholy[42])) return false;
	hodnoty[vrcholy[42]] = true;
	return true;
}

// Setting vertices that can be evaluate using base or complement of base, 12 vertices

inline bool Q6::Dopocitej_pom_baze()
{
  vrcholy[32] = magic_constant - vrcholy[16] - vrcholy[8] - vrcholy[4] - vrcholy[2] - vrcholy[1];
	if (pouzite(vrcholy[32]) || mimo_interval(vrcholy[32])) return false;
	hodnoty[vrcholy[32]] = true;

	vrcholy[33] = magic_constant - vrcholy[17] - vrcholy[9] - vrcholy[5] - vrcholy[3] - vrcholy[0];
	if (pouzite(vrcholy[33]) || mimo_interval(vrcholy[33])) return false;
	hodnoty[vrcholy[33]] = true;

	vrcholy[34] = magic_constant - vrcholy[18] - vrcholy[10] - vrcholy[6] - vrcholy[0] - vrcholy[3];
	if (pouzite(vrcholy[34]) || mimo_interval(vrcholy[34])) return false;
	hodnoty[vrcholy[34]] = true;

	vrcholy[35] = magic_constant - vrcholy[19] - vrcholy[11] - vrcholy[7] - vrcholy[1] - vrcholy[2];
	if (pouzite(vrcholy[35]) || mimo_interval(vrcholy[35])) return false;
	hodnoty[vrcholy[35]] = true;

	vrcholy[36] = magic_constant - vrcholy[20] - vrcholy[12] - vrcholy[0] - vrcholy[6] - vrcholy[5];
	if (pouzite(vrcholy[36]) || mimo_interval(vrcholy[36])) return false;
	hodnoty[vrcholy[36]] = true;

	vrcholy[37] = magic_constant - vrcholy[21] - vrcholy[13] - vrcholy[1] - vrcholy[7] - vrcholy[4];
	if (pouzite(vrcholy[37]) || mimo_interval(vrcholy[37])) return false;
	hodnoty[vrcholy[37]] = true;

	vrcholy[31] = 63 - vrcholy[32];
	if (pouzite(vrcholy[31])) return false;
	hodnoty[vrcholy[31]] = true;

	vrcholy[30] = 63 - vrcholy[33];
	if (pouzite(vrcholy[30])) return false;
	hodnoty[vrcholy[30]] = true;

	vrcholy[29] = 63 - vrcholy[34];
	if (pouzite(vrcholy[29])) return false;
	hodnoty[vrcholy[29]] = true;

	vrcholy[28] = 63 - vrcholy[35];
	if (pouzite(vrcholy[28])) return false;
	hodnoty[vrcholy[28]] = true;

	vrcholy[27] = 63 - vrcholy[36];
	if (pouzite(vrcholy[27])) return false;
	hodnoty[vrcholy[27]] = true;

	vrcholy[26] = 63 - vrcholy[37];
	if (pouzite(vrcholy[26])) return false;
	hodnoty[vrcholy[26]] = true;

	return true;
}

// Setting last vertices, 12 vertices

inline bool Q6::Dopocitej_zbyle()
{
  vrcholy[15] = magic_constant - vrcholy[46] - vrcholy[30] - vrcholy[6] - vrcholy[10] - vrcholy[12];
	if (pouzite(vrcholy[15]) || mimo_interval(vrcholy[15])) return false;
	hodnoty[vrcholy[15]] = true;

	vrcholy[14] = magic_constant - vrcholy[47] - vrcholy[31] - vrcholy[7] - vrcholy[11] - vrcholy[13];
	if (pouzite(vrcholy[14]) || mimo_interval(vrcholy[14])) return false;
	hodnoty[vrcholy[14]] = true;

	vrcholy[23] = magic_constant - vrcholy[54] - vrcholy[6] - vrcholy[30] - vrcholy[18] - vrcholy[20];
	if (pouzite(vrcholy[23]) || mimo_interval(vrcholy[23])) return false;
	hodnoty[vrcholy[23]] = true;

	vrcholy[22] = magic_constant - vrcholy[55] - vrcholy[7] - vrcholy[31] - vrcholy[19] - vrcholy[21];
	if (pouzite(vrcholy[22]) || mimo_interval(vrcholy[22])) return false;
	hodnoty[vrcholy[22]] = true;

	vrcholy[25] = magic_constant - vrcholy[56] - vrcholy[8] - vrcholy[16] - vrcholy[28] - vrcholy[26];
	if (pouzite(vrcholy[25]) || mimo_interval(vrcholy[25])) return false;
	hodnoty[vrcholy[25]] = true;

	vrcholy[24] = magic_constant - vrcholy[57] - vrcholy[9] - vrcholy[17] - vrcholy[29] - vrcholy[27];
	if (pouzite(vrcholy[24]) || mimo_interval(vrcholy[24])) return false;
	hodnoty[vrcholy[24]] = true;

	vrcholy[38] = 63 - vrcholy[25];
	if (pouzite(vrcholy[38])) return false;
	hodnoty[vrcholy[38]] = true;

	vrcholy[39] = 63 - vrcholy[24];
	if (pouzite(vrcholy[39])) return false;
	hodnoty[vrcholy[39]] = true;

	vrcholy[40] = 63 - vrcholy[23];
	if (pouzite(vrcholy[40])) return false;
	hodnoty[vrcholy[40]] = true;

	vrcholy[41] = 63 - vrcholy[22];
	if (pouzite(vrcholy[41])) return false;
	hodnoty[vrcholy[41]] = true;

	vrcholy[48] = 63 - vrcholy[15];
	if (pouzite(vrcholy[48])) return false;
	hodnoty[vrcholy[48]] = true;

	vrcholy[49] = 63 - vrcholy[14];
	return !pouzite(vrcholy[49]); // už nemusím nastavovat, že jsem použil hodnotu pro vrchol 49
}

// inner counting function - not exported

bool Q6::simuluj()
{
  if (!Nastav()) return false;
	if (!Dopocitej_kompl_baze()) return false;
	if (!Dopocitej_pom_baze()) return false;
	if (!Dopocitej_zbyle()) return false;
	return true;
}

// Exported counting function

// [[Rcpp::export]]
Rcpp::NumericVector count_dml()
{
  time_t t;
  srand(time(&t));
  init_gen();
  long long dml_counter = 0;
  long long iter_counter = 0;
  Q6 cube;
  cube.init();
  while (true)
  {
    if (cube.simuluj())
    {
      dml_counter++;
      cube.save_dml();
    }
    iter_counter++;
    try 
    {
      Rcpp::checkUserInterrupt();
    }
    catch (Rcpp::internal::InterruptedException e)
    {
      break;  
    }
  }
  NumericVector out(2);
  out[0] = dml_counter;
  out[1] = iter_counter;
  return out;
}