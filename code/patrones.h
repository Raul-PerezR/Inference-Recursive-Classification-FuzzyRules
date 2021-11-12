#ifndef _PATRONESH_
#define _PATRONESH_

#include <unordered_map>
#include <vector>
#include <string>
#include <list>
#include "VectorVar.h"
#include "example_set.h"
#include "vectordouble.h"


using namespace std;

struct TestResult{
  vector<double> cubiertos;
  vector<double> no_cubiertos;
  vector<double> acc;
  vector<double> error_intrinseco_porClase;
  double acierto_global;
  double acierto_sinNoCubiertos;
  double porcentaje_nuevos_patrones;
  double error_intrinseco;
};

void ProcesarResultados(TestResult & result);
void PintaResultadosTraining(const TestResult & result);
void PintaResultadosTest(const TestResult & result, bool conClases);
void CalcularMediaTestResult(vector<TestResult> lista, TestResult & result);
void SavePattern(const char *nomfich);

struct info{
  vector<double> conseq;
  vector<int> eje_list;
};

struct info2{
	double adaptacion;
	string subCadena;
};


struct infoUp {
	double key;
	int variable;
	int posCadenainicio;
	vector<info2> trozoAcambiar;
};



class Pattern{
private:
  unordered_map<string,info > diccionario;
	//diccionario.max_load_factor(3.0);
  int n_ejemplos;
  int n_clases;

  int BetterEtiqueta(const VectorVar &V, int variable, const example_set & E, int ejemplo);

public:
  Pattern(){
    diccionario.max_load_factor(3.0);
    n_ejemplos = 0;
    n_clases = 2;
  }

  Pattern( const VectorVar & V){
    diccionario.max_load_factor(3.0);
    n_ejemplos = 0;
    n_clases = V.SizeDomain(V.Consecuente());
  }

  Pattern (const Pattern & x);
  Pattern &operator=(const Pattern &x);



  void ExtraerPatronesBasicos(const example_set &Es, const VectorVar &V, TestResult & result);
  void ExtraerPatronesBasicosOriginalWM(const example_set &Es, const VectorVar &V, TestResult & result);
  void ExtraerPatronesBasicosAproximacionTFMRuben(const example_set &Es, const VectorVar &V, TestResult & result, int tries);

  void Aprendizaje_RecursivoUnEjemplo_WM_TFM_Ruben(const example_set &E, const VectorVar &V, const int eje,
  		 string cadena, int actualVar, double adapt, const vector<infoUp> &trozos, int &current_tries, int number_tries);


  void TestearPatronesBasicos(const example_set &Es, const VectorVar &V, TestResult & result);
  int  TestearUnEjemploPatronesBasicos(const example_set &Es, int eje, const VectorVar &V, double &peso);

  int InferenciaRecursivaOptimalizadaConProfundidadLimitada(const example_set &E, const VectorVar &V, const int eje, bool PCF, double &mu, string &antecedenteSeleccionado, int prof);
  void TestearRecursivoUnEjemploEficenteConProfundidadLimitada(const example_set &E, const VectorVar &V, const int eje,
  		 string cadena, int actualVar, double adapt, double &umbral, int &clase, bool PCF, const vector<infoUp> &trozos, string &antecedenteSeleccionado, int prof);

  int InferenciaRecursivaOptimalizadaConProfundidadLimitada_SalidaPorProfundidades(const example_set &E, const VectorVar &V, const int eje, bool PCF, double &mu, string &antecedenteSeleccionado, int prof);
  void TestearRecursivoUnEjemploEficenteConProfundidadLimitada_SalidaPorProfundidades(const example_set &E, const VectorVar &V, const int eje,
 		 string cadena, int actualVar, double adapt, vector<double> &umbral, vector<int> &clase, bool PCF, const vector<infoUp> &trozos, string &antecedenteSeleccionado,
 		 int prof, vector<pair<string,double> > &estructura);

  int InferenciaRecursivaOptimalizada(const example_set &E, const VectorVar &V, const int eje, bool PCF, double &mu, string &antecedenteSeleccionado, int TriesNumber);
  void TestearRecursivoUnEjemploEficente(const example_set &E, const VectorVar &V, const int eje,
		 string cadena, int actualVar, double adapt, double &umbral, int &clase, bool PCF, const vector<infoUp> &trozos, string &antecedenteSeleccionado, int &TriesNunmber);

  void TestearRecursivoUnEjemplo(const example_set &E, const VectorVar &V, const int eje, string cadena, int actualVar, double adapt, double &umbral, int &clase, bool PCF, string &antecedenteSeleccionado);
  void TestearRecursivo(const example_set &Es, const VectorVar &V, TestResult & result, bool PCF, int TriesNumber);
  void TestearRecursivoUnEjemplo(const example_set &E, const VectorVar &V, const int eje, string cadena, int actualVar, double adapt, double &umbral, int &clase, bool PCF, string &antecedenteSeleccionado, int &TriesNumber);

  int TesteoDistanciaHamming(const example_set &Es, const VectorVar &V, const int eje);
  void TestearRecursivoVariosDicionarios(const example_set &Es, const VectorVar &V, const Pattern &P3, const VectorVar &V3, const Pattern &P2, const VectorVar &V2, TestResult & result);

  int TestearPatronesBasicosClassicOriginal_UnEjemplo(const example_set &Es, const VectorVar &V, const int eje, const vector<double> &pesoR, const vector<int> &claseR, double &mu, string &antecedente);
  int TestearPatronesBasicosClassic_UnEjemplo(const example_set &Es, const VectorVar &V, const int eje, const vector<double> &pesoR, const vector<int> &claseR, double &mu, string &antecedente);
  void TestearPatronesBasicosClassic(const example_set &Es, const VectorVar &V, TestResult & result);
  void TestearPatronesBasicosClassicDisparos(const example_set &Es, const VectorVar &V, TestResult & result, vector<pair<int, pair <string, double> > > &disparos);
  void CalcularPesoYClases(vector<double> &pesoR, vector<int> &claseR,  bool PCF);



  int N_Pattern(){return diccionario.size();}
  int N_Patrones_paraClase (int clase) const;
  int N_Ejemplos_paraClase (int clase) const;

  int InferirDiccionarioPatrones (vectordouble &ejemplo, const VectorVar &V);
  double InferirDicionario (const example_set &E, const VectorVar &V);


};



#endif
