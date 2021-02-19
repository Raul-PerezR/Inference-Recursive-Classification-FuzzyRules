#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/times.h>
#include <math.h>
#include <fstream>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <queue>
#include <omp.h>
#include "domain_t.h"
#include "VectorVar.h"
#include "vectordouble.h"
#include "example_set.h"
#include "multipoblacion.h"
#include "genetcode.h"
#include "ruleset.h"
#include "Relaciones.h"
#include "operaciones.h"
#include "patrones.h"

using namespace std;


double Weight_Individual(multipoblacion &G, const VectorVar &V, const example_set &E, VectorVar &F, double ***Tabla, double ***Tabla_Nuevas, int elemento_i,
			 double &n_mas, double &n_menos);

void ReajustarPesoaReglas(multipoblacion &G, const VectorVar &V, example_set &E, ruleset &R, VectorVar &F,
			  double ***Tabla, double ***Tabla_Nuevas);

int GA(multipoblacion &G, const VectorVar &V, example_set &E, ruleset &R, ruleset &R_dis, ruleset &R_mix, double **I, double *agregado,
       double *peso_agregado, int Total_ejemplos, bool discriminante, bool primeraVez, bool con_peso,
       int *ponderacion_por_clase, double ***Tabla, double ***Tabla_Nuevas, VectorVar &F, int *basica_para_regla,
       int *errores_para_regla, int *reglas_a_eliminar, int &n_reglas_a_eliminar, bool &Added_Rule, bool &Specific_Rule,
       int tam_rule_not_removed);

//void mGA(multipoblacion *G, int tam, const VectorVar &V, example_set &E);
void fitness (multipoblacion &G, const VectorVar &V, const example_set &E, ruleset &R, bool discriminante, bool con_peso,
              int *ponderacion_por_clase, int nrule, double ***Tabla, double ***Tabla_Nuevas, VectorVar &F, int *basica_para_regla,
	      int *errores_para_regla, int tam_rule_not_removed);

void Obtener_Regla_Discriminante(multipoblacion &G, const VectorVar &V, example_set &E, ruleset &R, int elemento_i,
                               bool discriminante, bool reducir_reglas, int* ponderacion_por_clase, double ***Tabla,
                               double ***Tabla_Nuevas, VectorVar &F, int *basica_para_regla,
			       int *errores_para_regla, int *reglas_a_eliminar, int &n_reglas_a_eliminar, bool debug,
			       string patr, bool inGA, string initial_rule, int tam_rule_not_removed);


void Obtener_Regla_Discriminante_E(multipoblacion &G, const VectorVar &V, example_set &E, ruleset &R, int elemento_i,
                               bool discriminante, bool reducir_reglas, int* ponderacion_por_clase, double ***Tabla,
                               double ***Tabla_Nuevas, VectorVar &F, int *basica_para_regla,
			       int *errores_para_regla, int *reglas_a_eliminar, int &n_reglas_a_eliminar, bool debug, string patr,
			       bool inGA, string initial_rule, double nmas, double nmenos);

string Patron1Regla(ruleset &R_1, const VectorVar &V, int pos);

string Patron2Patrones(string pat1, string pat2);

bool ComprobacionEspGen(string esp, string gen);


string Patron2Reglas(ruleset &R_1, ruleset &R_2, const VectorVar &V, int pos);

void Pinta2Reglas(ruleset &R_1, ruleset &R_2, const VectorVar &V, int pos);



void Obtener_Regla_Descriptiva(multipoblacion &G, const VectorVar &V, example_set &E, ruleset &R, int elemento_i,
                               bool discriminante, bool reducir_reglas, int* ponderacion_por_clase, double ***Tabla,
                               double ***Tabla_Nuevas, VectorVar &F, int *basica_para_regla,
			       int *errores_para_regla, int *reglas_a_eliminar, int &n_reglas_a_eliminar, bool debug,
			       string patr, bool inGA, string initial_rule, int tam_rule_not_removed);


void Obtener_Regla_Descriptiva_E(multipoblacion &G, const VectorVar &V, example_set &E, ruleset &R, int elemento_i,
                               bool discriminante, bool reducir_reglas, int* ponderacion_por_clase, double ***Tabla,
                               double ***Tabla_Nuevas, VectorVar &F, int *basica_para_regla,
			       int *errores_para_regla, int *reglas_a_eliminar, int &n_reglas_a_eliminar, bool debug, string patr,
			       bool inGA, string initial_rule, double nmas, double nmenos);


double Bondad (multipoblacion &G, const VectorVar &V, const example_set &E, ruleset &R, double *valor, int N_valor,
               int elemento_i, bool discriminante, bool con_peso, int *ponderacion_por_clase, int nrule,
               double ***Tabla, double ***Tabla_Nuevas, VectorVar &F, int *basica_para_regla, int *errores_para_regla, int *reglas_a_eliminar, int &n_reglas_a_eliminar,
	       bool esp_gen, string initial_rule, int tam_rule_not_removed);

double Bondad_E (multipoblacion &G, const VectorVar &V, const example_set &E, ruleset &R, double *valor, int N_valor,
               int elemento_i, bool discriminante, bool con_peso, int *ponderacion_por_clase, int nrule,
               double ***Tabla, double ***Tabla_Nuevas, VectorVar &F, int *basica_para_regla, int *errores_para_regla, int *reglas_a_eliminar, int &n_reglas_a_eliminar,
	       bool esp_gen, string initial_rule, double nmas, double nmenos);


double Bondad_Final (multipoblacion &G, const VectorVar &V, const example_set &E, ruleset &R, double *valor, int N_valor,
               int elemento_i, bool discriminante, bool con_peso, int *ponderacion_por_clase, int nrule,
               double ***Tabla, double ***Tabla_Nuevas, VectorVar &F, int *basica_para_regla, int *errores_para_regla, int *reglas_a_eliminar, int &n_reglas_a_eliminar);

example_set RetenerEjemplos(ruleset &R, const VectorVar &V, example_set &E, VectorVar &F);


double Aciertos(const ruleset &R, const VectorVar &V, example_set &E, double *agregado, double *peso_agregado,
                string nomfich, VectorVar &F);

double Aciertos(ruleset &R, const VectorVar &V, example_set &E, VectorVar &F);
double AciertosLigero(ruleset &R, const VectorVar &V, example_set &E, VectorVar &F);

double Aciertos(ruleset &R, const VectorVar &V, example_set &E, string nomfich, VectorVar &F);

double Filtrar_Reglas(ruleset &R, const VectorVar &V, example_set &E, VectorVar &F);
double Filtrar_Reglas(ruleset &R, const VectorVar &V, example_set &E, ruleset &R_Descrip, VectorVar &F);
double Filtrar_Reglas_New(ruleset &R, ruleset &R_dis, ruleset &R_mix, int *reglas_a_eliminar,
			  int n_reglas_a_eliminar, bool Specific_Rule);



double Inferencia2BasesReglas(ruleset &R, ruleset &R_Descrip, const VectorVar &V, example_set &E, int *fallosR, int *fallosR_Descrip, int *disparosR, int *disparosR_Descrip);

void Error(const ruleset &R, const example_set &E);

void PintaIndividuo(genetcode GC, const VectorVar &V);

double **ReservaParaMedidasDeInformacion(const VectorVar &V);
void LiberaMedidasDeInformacion(double **&matriz, int n);
void MedidasDeInformacion(const VectorVar &V, const example_set &E, double **I);

double CalcularInformacionCondicionada2variablesRespectoClase (const VectorVar &V, const example_set &E,
        int casos1, int casos2, int casosclase,
        double ***Conjunta123, double **Conjunta13,
        double **Conjunta23, double *simple3, double *I);


void CalcularProbabilidadConjunta3Variables(const VectorVar &V, const example_set &E,
        int var1, int var2, int var3, int n_casos1,
        int n_casos2, int n_casos3, double ***m);

string extract_noun(string file);

bool Existe(string nom_fich);

void Max_Min(const example_set &E, const VectorVar &V, int v1, int v2, int op, double &min, double &max);

void MedidasDeInformacion_operadores(VectorVar &V, example_set &E, VectorVar &F);
void MedidasDeInformacion_operadores_AjustandoDominio(VectorVar &V, example_set &E, VectorVar &F);

void Redefiniendo_Particiones_Variables (VectorVar &V, example_set &E);
void MedidasDeInformacion_operandoresUnarios(VectorVar &V, example_set &E, int nuevo,bool *bina);
void MedidasDeInformacion_operandores2(VectorVar &V, example_set &E, int nuevos, bool *bina);
double CalcularInformacion (int tamaX, int tamaY, double *px, double *py, double **pxy);
double CalcularInformacion (int tamaX, int tamaY, double *px, double *py, double **pxy, double *I);
string Simbolo_Operador(int op);
void MedidasDeInformacion_2variables_cond_consecuente(VectorVar &V, example_set &E);

bool Terminado_Aprendizaje(multipoblacion &G, example_set &E, bool &ultimo_modo);
void Obtener_Aciertos_y_Errores_Reglas_Aprendidas(example_set &E, int *basica_para_regla, int *errores_para_regla, int nrule);

double Aciertos_Lista_Reglas_Eliminables(ruleset &R, const VectorVar &V, example_set &E, VectorVar &F, int*lista, int &tama, int *l_aciertos, int &tama_l);




//---- Funciones para manejar la tabla con las adaptaciones de los ejemplos a las etiquetas de las variables
double ***Crear_Tabla_Adaptaciones(const VectorVar &V, const example_set &E);
double ***Crear_Tabla_Adaptaciones_Nuevas(const VectorVar &V, const example_set &E);
void Destruir_Tabla_Adaptaciones(const VectorVar &V, const example_set &E, double ***&Tabla);
double T_Adaptacion_Variable_Etiqueta(const VectorVar &V, double ***Tabla, int ejemplo,
                                      int variable, int etiqueta);
double T_Adaptacion_Variable(const VectorVar &V, double ***Tabla, int ejemplo, int variable, string labels);
double T_Adaptacion_Antecedente(const VectorVar &V, double ***Tabla, int ejemplo, string regla,
                                double *nivel_variable, double umbral);
double T_Adaptacion_Antecedente(const VectorVar &V, double ***Tabla, int ejemplo, string regla, double umbral);
double T_Adaptacion_Antecedente_Nuevas(const VectorVar &V, double ***Tabla, int ejemplo, string regla,
                                       int *nivel_variable_funciones, int tama_nivel_variable_funciones, double matching_previo);
void T_AdaptacionC(const VectorVar &V, double ***Tabla, int ejemplo, int etiq, double &pos, double &neg);
int BetterEtiqueta(const VectorVar &V,double ***Tabla, int variable, int ejemplo);



double Evento1_Especificar(multipoblacion &G, const VectorVar &V, example_set &E_Par, ruleset &R, ruleset &R_dis, ruleset &R_mix, int *ponderacion_por_clase, VectorVar &F,
       double ***Tabla, double ***Tabla_Nuevas, int regla, double *multi_valor);

double Evento1_Generalizar(multipoblacion &G, const VectorVar &V, example_set &E_Par, ruleset &R, ruleset &R_dis, ruleset &R_mix, int *ponderacion_por_clase, VectorVar &F,
       double ***Tabla, double ***Tabla_Nuevas, int regla, double *multi_valor);

void PintaDiferenciaPatrones(const string &cad1, const string &cad2);
int Comparar2patrones(const string &cad1, const string &cad2);
int Comparar2patronesPinta(const string &cad1, const string &cad2);
bool AnalizarConjuntosDeReglas(const VectorVar &V, ruleset &R, ruleset &R_Descrip, ruleset &R_mix, int nuevas);
void UnionReglasDescriptivas(ruleset &R1, int r1, ruleset &R2, int r2);
void IntersecionReglasDiscriminantes(ruleset &R1, int r1, ruleset &R2, int r2);

string Interseccion(const string &cad1, const string &cad2);
string Union(const string &cad1, const string &cad2);
int Grados_Libertad(string cad);
int DistanciaHamming(const string &s1, const string &s2);
int DiferenciaDistanciaHamming(const string &s1, const string &s2, int pos_inicial);


void ActualizarPeso(multipoblacion &G, const VectorVar &V, example_set &E, ruleset &R, VectorVar &F,
       double ***Tabla, double ***Tabla_Nuevas);


bool vectoresIguales_E(double *v1, double *v2, int tama);
bool vectoresIguales(double *v1, double *v2, int tama);

void Pausa(){
  cout << "Pulsa una Tecla para continuar...\n";
  char ch;
  cin >> ch;
}

string itos(int n)
{
   const int max_size = 50;
   char buffer[max_size] = {0};
   sprintf(buffer, "%d", n);
   return string(buffer);
}

string GeneraString(int pos, int tama){
	string aux="";
	for (int i=0; i<tama; i++)
	  aux.push_back('0');
	aux[pos]='1';
	return aux;
}

vector<char> simbolo;

typedef pair<double, double> pareja;

struct ComparaMultiMap{
	bool operator()(const pareja &s1, const pareja &s2){
		return (s1.first<s2.first or
			      (s1.first==s2.first and s1.second<s2.second)
					 );
	}
};


string InterpretaRegla(const VectorVar &V, string regla){
	// Interpretar la regla
	string jugada;
	int longitud=0;
	for (int j=0; j<V.N_Antecedente(); j++){
			if (regla[longitud]=='1' and regla[longitud+1]=='1'){
				jugada.push_back('-');
			}
			else {
				int l=0;
	  		while (regla[longitud+l]!='1')
					l++;
  			jugada.push_back(simbolo[l]);
			}
			longitud+=V.SizeDomain(j);
	}
	return jugada;
}

void SalvarMap(const unordered_map<string, vector<int> > & Conjunto){
	ofstream f ("MapFile.dat");
	for (auto it=Conjunto.begin(); it!=Conjunto.end(); it++){
		f << it->first;
		f << " " << it->second.size();
		for (int i=0; i<it->second.size(); i++)
		   f << " " << it->second[i];
		f << endl;
	}
	f.close();
}

unordered_map<string, vector<int> > CargarMap(){
	unordered_map<string, vector<int> > Conjunto;
	Conjunto.max_load_factor(3.0);

	ifstream f("MapFile.dat");
	string clave;
	int tama, aux;
	vector<int> v;
	while (!f.eof()){
		v.clear();
		f >> clave;
		f >> tama;
		for (int i=0; i<tama; i++){
			f >> aux;
			v.push_back(aux);
		}

		Conjunto.insert(pair<string,vector<int> > (clave,v));
	}

	return Conjunto;
}


void SalvarSet(const set<string> & Conjunto){
	ofstream f ("SetFile.dat");
	for (auto it=Conjunto.begin(); it!=Conjunto.end(); it++)
		f << *it << endl;
	f.close();
}

set<string> CargarSet(){
	set<string> Conjunto;

	ifstream f("SetFile.dat");
	string clave;
	int tama, aux;
	while (!f.eof()){
		f >> clave;
		Conjunto.insert(clave);
	}

	return Conjunto;
}





string Componer(const VectorVar &V, string regla, string aux, int var){
	int i=0, longitud=0;
	while (i<var){
		longitud+=V.SizeDomain(i);
		i++;
	}
	for (int l=0; l<V.SizeDomain(var);l++)
		regla[longitud+l]=aux[l];
	return regla;
}




bool Incluido(string s1, string s2){
	string s3=s1;
	for (int i=0; i<s1.length(); i++){
		if (s1[i]=='1' and s2[i]=='1')
			s3[i]='1';
		else
		  s3[i]='0';
  }

	/*if (s1==s3){
	cout << "Comparando: " << s1 << endl;
	cout << "Comparando: " << s2 << endl;
	cout << "          : " << s3 << endl << endl << endl;
}*/

	return (s3==s1);
}

string InterseccionN(string s1, string s2){
	string s3=s1;
	for (int i=0; i<s1.length(); i++){
		if (s1[i]=='1' and s2[i]=='1')
			s3[i]='1';
		else
		  s3[i]='0';
  }
  return s3;
}


bool IncluidoOnSet(const set<string> & Conjunto, string s){
	bool on=false;

  for (auto it=Conjunto.begin(); it!=Conjunto.end() and !on; it++){
		on = Incluido(s, *it);
	}

  return on;
}

bool IncluidoOnSet(const unordered_map<string, vector<int> > & Conjunto, string s){
	bool on=false;

  for (auto it=Conjunto.begin(); it!=Conjunto.end() and !on; it=Conjunto.end())
		on = Incluido(s, it->first);

  return on;
}



bool Infer(const example_set &Es, const VectorVar &V, double ***Tabla,
				   string regla, double *varlevel, int &clase0, int &clase1, bool &salir){

	double adaptC1;
	clase0=clase1=0;
	double *adaptA = new double[Es.N_Examples()];
	double *adaptC0 = new double[Es.N_Examples()];


  //#pragma omp parallel for num_threads(4)
	for (int e=0; e<Es.N_Examples() and !salir; e++){
		if (!salir){
			adaptA[e] = T_Adaptacion_Antecedente(V,Tabla,e,regla,varlevel,0.5);
			T_AdaptacionC(V,Tabla,e,0,adaptC0[e],adaptC1);
			if (adaptA[e]>0 and adaptC0[e]>0){
				#pragma omp atomic
				clase0++;
			}
			else if (adaptA[e]>0){
				#pragma omp atomic
				clase1++;
			}
		}
	}

  delete [] adaptA;
	delete [] adaptC0;

	if (salir)
	  return false;
	else
	  return true;
	//varlevel[j]=1;
	//cout << j << ") clase 0: " << clase0 << "  clase 1: " << clase1 << endl;
}


bool CubreEjemplos(string regla,const VectorVar &V, const example_set &Es, double ***Tabla, vector<bool> & cubierto){
	double *varlevel = new double[V.N_Antecedente()+1];
	for (int i=0; i< V.N_Antecedente(); i++){
			varlevel[i]=1;
	}
	varlevel[V.N_Antecedente()]=0.5;

  bool modificado=false;

	double adaptA;

	for (int e=0; e<Es.N_Examples(); e++){
		if (!cubierto[e]){
			adaptA = T_Adaptacion_Antecedente(V,Tabla,e,regla,varlevel,0.5);
			if (adaptA>0){
				modificado = true;
				cubierto[e]=true;
			}
		}
	}

	delete [] varlevel;

	return modificado;
}



bool CalSobreMap(const VectorVar &V, string regla, const unordered_map<string, vector<int> > & dic ,
	const unordered_map<string, vector<int> > & dic2 ,int &clase0, int &clase1,
  const example_set &E, double ***Tabla, int var, bool &salir){
	set<string> Conjunto, ConjuntoVar, ConjuntoAcum, ConjuntoComp;
	string *valores = new string[V.N_Antecedente()];
	int i=0;
	//cout << "Regla: " << InterpretaRegla(V,regla) << endl;
	for (int j=0; j<V.N_Antecedente(); j++){
		valores[j]=regla.substr(i+0,V.SizeDomain(j));
		//cout << "Variable: " << valores[j] << endl;
		i+=V.SizeDomain(j);
	}

  ConjuntoAcum.insert("");

	for (int i=0; i<V.N_Antecedente() and !salir; i++){

		for (int l=0; l<V.SizeDomain(i) and !salir; l++){
		  if (valores[i][l]=='1'){
				string aux = GeneraString(l,V.SizeDomain(i));
				string nrule = Componer(V, regla, aux, i);
				auto it=dic2.find(nrule);
				if (it!=dic2.end()){
				  //cout << "\nnregla: " << nrule << endl;
					ConjuntoComp.insert(nrule);
				}
				else{
				  //cout << aux << endl;
				  ConjuntoVar.insert(aux);
				}
			}
		}

		for (auto it2=ConjuntoVar.begin(); it2!=ConjuntoVar.end() and !salir; it2++){
			for (auto it=ConjuntoAcum.begin(); it!=ConjuntoAcum.end() and !salir; it++){
				string aux=*it+*it2;
				//cout << aux << endl;
				Conjunto.insert(aux);
			}
		}

		ConjuntoVar.clear();
		if (i!=V.N_Antecedente()-1){
			ConjuntoAcum = Conjunto;
			Conjunto.clear();
		}
	}

	// Calculo

	//cout << "Tama Reglas compues: " << ConjuntoComp.size() << endl;

  int t=0;
	clase0=0;
	clase1=0;

	for (auto it=ConjuntoComp.begin(); it!=ConjuntoComp.end() and !salir; it++){
		auto it2 = dic2.find(*it);
		if (it2!=dic2.end()){
			//cout << "------ " << InterpretaRegla(V,*it) << "\t" << it2->second[0] << "\t" << it2->second[1] << endl;
		  clase0+=it2->second[0];
			clase1+=it2->second[1];
			t++;
		}
	}

	for (auto it=ConjuntoComp.begin(); it!=ConjuntoComp.end() and !salir; it++){
		auto it2 = dic2.find(*it);
		for (auto it3= ++it; it3!=ConjuntoComp.end() and !salir; it3++){
			auto it4 = dic2.find(*it3);
			if (it4!=it2 and it4!=dic2.end()){
				string inter = InterseccionN(it2->first, it4->first);
				int c0=0, c1=0;
				auto it5 = dic.find(inter);
				if (it5!=dic.end()){
					c0 = it5->second[0];
					c1 = it5->second[1];
				}
				else {
					auto it6 = dic2.find(inter);
					if (it6!=dic2.end()){
						c0 = it6->second[0];
						c1 = it6->second[1];
					}
					else {
						//cout << inter << " hay que calcularlo\n";
						salir = true;
					}
				}
				//CalSobreMap(V, inter, dic, dic2 ,c0,c1,E,Tabla,var);
				//cout << "InterS " << InterpretaRegla(V,inter) << "\t" << c0 << "\t" << c1 << endl;
				clase0-=c0;
				clase1-=c1;
			}
		}
	}


  //cout << "Tama Reglas simples: " << Conjunto.size() << endl;


	for (auto it=Conjunto.begin(); it!=Conjunto.end() and !salir; it++){
		if (!IncluidoOnSet(ConjuntoComp, *it)){
			auto it2 = dic.find(*it);
			if (it2!=dic.end()){
				//cout << "       " << InterpretaRegla(V,*it) << "\t" << it2->second[0] << "\t" << it2->second[1] << endl;
				clase0+=it2->second[0];
				clase1+=it2->second[1];
			}
			else {
				//cout << "       " << InterpretaRegla(V,*it) << "\t" << 0 << "\t" << 0 << endl;
			}
		}
	}

	//cout << "Suma de todo es    : " << clase0 << "\t" << clase1 << endl;

  /*double *vlevel = new double[V.N_Antecedente()+1];
  for (int j=0; j<V.N_Antecedente(); j++)
	  vlevel[j]=1;
	vlevel[V.N_Antecedente()]=0.5;

	int nc0, nc1;
	Infer(E,V,Tabla,regla,vlevel,nc0,nc1);

	//cout << "\nEl correcto es     : " << nc0 << "\t" << nc1 << endl;
	//cout <<   "Obtenido           : " << clase0 << "\t" << clase1 << endl;

	delete [] vlevel;

  if (nc0!=clase0 or nc1!=clase1)
  	Pausa();
  */
  //cout << endl << endl;
	delete [] valores;

	return !salir;

}


void deleteLine(){
	putchar(27);    //codigo ascii decimal de ESC
	putchar('[');
	putchar('0');    //cero para borrar desde el cursor hasta el final
	putchar('K');

	putchar(27);    //codigo ascii decimal de ESC
	putchar('[');
	putchar('1');   // Se mueve una linea arriba
	putchar('A');

	putchar(27);    //codigo ascii decimal de ESC
	putchar('[');
	putchar('0');    //cero para borrar desde el cursor hasta el final
	putchar('K');
}



//-----------------------------------------------------------------------------------------------------


void OrdenarPatrones( const VectorVar &V,
	                    unordered_map<string,vector<int> > & diccionario,
				      			  multimap <pareja, unordered_map<string,vector<int> >::iterator, ComparaMultiMap > &reves){


	for (auto it=diccionario.begin(); it!=diccionario.end(); it++){
		pareja sum;
		sum.first=0;

		int mayor=0;
		for (int k=0; k<V.SizeDomain(V.Consecuente()); k++){
			if (it->second[k]>it->second[mayor])
			mayor=k;
		}

		for (int k=0; k<V.SizeDomain(V.Consecuente()); k++){
			//if (k!=mayor)
			 sum.first+= it->second[k];
		}
		//sum.first = it->second[mayor]/sum.first;
		sum.second = it->second[mayor];

		reves.insert(pair<pareja, unordered_map<string,vector<int> >::iterator >  (sum, it) );
	}

	//return reves;
}

//-----------------------------------------------------------------------------------------------------
double CalPeso(vector<int> &v, int &clase){
	double suma=0, peso;
	int mayor=0;
	for (int i=1; i<v.size(); i++){
		suma += v[i];
		if (v[i]>v[mayor]){
			mayor = i;
		}
	}
	if (suma>0){
		clase = mayor;
		peso = mayor / (suma);
	}
	else{
		clase = -1;
		peso = -1;
	}
}

//-----------------------------------------------------------------------------------------------------

	int  InferirDiccionarioPatrones (int ejemplo, const VectorVar &V, double ***Tabla, unordered_map<string,vector<int> > & diccionario){
		auto it = diccionario.begin();
		double adapt_aux, adaptacion = 0, peso;
		int clase = -1;

		while (it!=diccionario.end() and adaptacion < 1 ){
			adapt_aux = T_Adaptacion_Antecedente(V,Tabla,ejemplo,it->first,0);
			if (adapt_aux>0){
				peso = CalPeso(it->second, clase);
				if (adapt_aux > adaptacion){
					adaptacion = adapt_aux;
				}
			}
			it++;
		}
		return clase;
	}



//-----------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------

void  SacarPatronesBasicos(const example_set &Es, const VectorVar &V, double ***Tabla, unordered_map<string,vector<int> > & diccionario, int ejemplo){
	//Meter los patrones en un map
	//unordered_map<string,vector<int> > diccionario;
	//diccionario.max_load_factor(3.0);


	// Construyo el patron
	for (int i=0; i<Es.N_Examples(); i++){
		string aux="";
		for (int j=0; j<V.N_Antecedente(); j++){
			string auxaux="";
			for (int l=0; l<V.SizeDomain(j); l++){
				auxaux.push_back('0');
			}
			//auxaux.push_back(' ');
			auxaux[BetterEtiqueta(V,Tabla,j,i)]='1';
			aux += auxaux;
		}

		auto it = diccionario.find(aux);
		if (it==diccionario.end()){
			double adaptacion, peso;
			int clase;
			int claseDelEjemplo = BetterEtiqueta(V, Tabla, V.Consecuente(), i);
			clase = InferirDiccionarioPatrones(i, V, Tabla, diccionario);
			if (claseDelEjemplo!=clase){
				vector<int> v;
				for (int k=0; k<V.SizeDomain(V.Consecuente())+1; k++)
				v.push_back(0);
				v[BetterEtiqueta(V,Tabla,V.Consecuente(),i)]++;
				diccionario.insert(pair<string,vector<int> > (aux,v));
			}
		}
		else {
			it->second[BetterEtiqueta(V,Tabla,V.Consecuente(),i)]++;
		}
	}


	// Datos de los patrones Encontrados

	//return diccionario;
}


//-----------------------------------------------------------------------------------------------------

  unordered_map<string,vector<int> >
  SacarPatronesBasicos(const example_set &Es, const VectorVar &V, double ***Tabla){
    //Meter los patrones en un map
    unordered_map<string,vector<int> > diccionario;
    diccionario.max_load_factor(3.0);

    for (int i=0; i<Es.N_Examples(); i++){
      string aux="";
      for (int j=0; j<V.N_Antecedente(); j++){
        string auxaux="";
        for (int l=0; l<V.SizeDomain(j); l++){
          auxaux.push_back('0');
        }
        //auxaux.push_back(' ');
        auxaux[BetterEtiqueta(V,Tabla,j,i)]='1';
        aux += auxaux;
      }

      auto it = diccionario.find(aux);
      if (it==diccionario.end()){
        vector<int> v;
        for (int k=0; k<V.SizeDomain(V.Consecuente())+1; k++)
          v.push_back(0);
        v[BetterEtiqueta(V,Tabla,V.Consecuente(),i)]++;
        diccionario.insert(pair<string,vector<int> > (aux,v));
      }
      else {
        it->second[BetterEtiqueta(V,Tabla,V.Consecuente(),i)]++;
      }
    }

    // Datos de los patrones Encontrados

    return diccionario;
  }

	//-----------------------------------------------------------------------------------------------------

	  void  SacarPatronesBasicos(const example_set &Es, const VectorVar &V, double ***Tabla, unordered_map<string,vector<int> > & diccionario){
	    //Meter los patrones en un map
	    //unordered_map<string,vector<int> > diccionario;
	    //diccionario.max_load_factor(3.0);


			// Construyo el patron
	    for (int i=0; i<Es.N_Examples(); i++){
	      string aux="";
	      for (int j=0; j<V.N_Antecedente(); j++){
	        string auxaux="";
	        for (int l=0; l<V.SizeDomain(j); l++){
	          auxaux.push_back('0');
	        }
	        //auxaux.push_back(' ');
	        auxaux[BetterEtiqueta(V,Tabla,j,i)]='1';
	        aux += auxaux;
	      }

	      auto it = diccionario.find(aux);
	      if (it==diccionario.end()){
	        vector<int> v;
	        for (int k=0; k<V.SizeDomain(V.Consecuente())+1; k++)
	          v.push_back(0);
	        v[BetterEtiqueta(V,Tabla,V.Consecuente(),i)]++;
	        diccionario.insert(pair<string,vector<int> > (aux,v));
	      }
	      else {
	        it->second[BetterEtiqueta(V,Tabla,V.Consecuente(),i)]++;
	      }
	    }

	    // Datos de los patrones Encontrados

	    //return diccionario;
	  }


	//-----------------------------------------------------------------------------------------------------

	int SacarPatronesBasicosTest(const example_set &Es, const VectorVar &V, double ***Tabla,  unordered_map<string,vector<int> > & diccionario, int &no_cubiertos){
	    //Meter los patrones en un map

      no_cubiertos=0;
	    diccionario.max_load_factor(3.0);

			int acc=0;

	    for (int i=0; i<Es.N_Examples(); i++){
	      string aux="";
	      for (int j=0; j<V.N_Antecedente(); j++){
	        string auxaux="";
	        for (int l=0; l<V.SizeDomain(j); l++){
	          auxaux.push_back('0');
	        }
	        //auxaux.push_back(' ');
	        auxaux[BetterEtiqueta(V,Tabla,j,i)]='1';
	        aux += auxaux;
	      }

	      auto it = diccionario.find(aux);
				int mayor;
	      if (it!=diccionario.end()){
					mayor=0;
		      for (int k=1; k<V.SizeDomain(V.Consecuente()); k++){
		        if (it->second[k]>=it->second[mayor])
		          mayor=k;
		      }
					if (mayor==Es.Data(i,V.Consecuente()))
	        	acc++;
	      }
				else{
						no_cubiertos++;
				}

	    }
	    // Datos de los patrones Encontrados
	    return acc;
	  }

//-----------------------------------------------------------------------------------------------------

		int SacarPatronesBasicosTest(const example_set &Es, const VectorVar &V, double ***Tabla,  unordered_map<string,vector<int> > & diccionario, int &no_cubiertos, int ejemplo){
		    //Meter los patrones en un map

	      no_cubiertos=0;
		    diccionario.max_load_factor(3.0);

				int acc=0;

		    for (int i=0; i<Es.N_Examples(); i++){
		      string aux="";
		      for (int j=0; j<V.N_Antecedente(); j++){
		        string auxaux="";
		        for (int l=0; l<V.SizeDomain(j); l++){
		          auxaux.push_back('0');
		        }
		        //auxaux.push_back(' ');
		        auxaux[BetterEtiqueta(V,Tabla,j,i)]='1';
		        aux += auxaux;
		      }

		      auto it = diccionario.find(aux);
					int mayor;
		      if (it!=diccionario.end()){
						mayor=0;
			      for (int k=1; k<V.SizeDomain(V.Consecuente()); k++){
			        if (it->second[k]>=it->second[mayor])
			          mayor=k;
			      }
						if (mayor==Es.Data(i,V.Consecuente()))
		        	acc++;
		      }
					else{
						// El patron no existe
						double adaptacion, peso;
						int clase;
						int claseDelEjemplo = BetterEtiqueta(V, Tabla, V.Consecuente(), i);
						clase = InferirDiccionarioPatrones(i, V, Tabla, diccionario);
						if (clase != -1){
							if (clase==Es.Data(i,V.Consecuente()))
			        	acc++;
						}
						else {
							no_cubiertos++;
						}
					}

		    }
		    // Datos de los patrones Encontrados
		    return acc;
		  }


//-----------------------------------------------------------------------------------------------------

  double CalculoReglas(const example_set &Es, const VectorVar &V, string regla2, double ***Tabla,
									   int j, double *varlevel,
                     const unordered_map<string, vector<int> > & diccionario,
                     unordered_map<string, vector<int> > & diccionario2,
                     int &clase0, int &clase1){
      int inCero, inUno;
      bool section1=false, section2=false;
			double t0, t1;

			t0 = omp_get_wtime();

      #pragma omp parallel num_threads(2) shared (V,regla2, Es, Tabla, section1, section2)
      {
        #pragma omp sections nowait
        {

          #pragma omp section
          {
            //int ID = omp_get_thread_num();
            //cout << "\t\tEmpieza Cal\n";
            section1 = CalSobreMap(V,regla2,diccionario,diccionario2,clase0,clase1, Es, Tabla, j, section2);
            //cout << "\t\tCal: " << clase0 << " " << clase1 << " Procesador: " << ID << endl;
            //section1=true;
          }

          #pragma omp section
          {
            //int ID = omp_get_thread_num();
            //cout << "\t\t\tEmpieza Inf\n";
            section2 = Infer(Es, V, Tabla, regla2, varlevel, inCero, inUno, section1);
            //cout << "\t\t\tInf: " << inCero << " " << inUno << " Procesador: " << ID << endl;
            //section2=true;
          }

        }
      }
			t1 = omp_get_wtime();
      if (section2){
        clase0 = inCero;
        clase1 = inUno;
      }

			//double tconsumido = 1.0*(t11-t10)/CLOCKS_PER_SEC;
			double tconsumido = t1-t0;

      vector<int> v;
      v.push_back(clase0);
      v.push_back(clase1);
      diccionario2.insert(pair<string,vector<int> > (regla2,v));

			return tconsumido;
  }


  //-----------------------------------------------------------------------------------------------------


   map<string, vector<double> >
   ReducirSet (const example_set &Es, const VectorVar &V, double ***Tabla,
              unordered_map<string,vector<int> > &diccionario,
              unordered_map<string,vector<int> > &diccionario2,
              set<string> &rFinal){

    multimap<pareja, unordered_map<string,vector<int> >::iterator, ComparaMultiMap > reves;

    for (auto it=rFinal.begin(); it!=rFinal.end(); it++){
      auto it2 = diccionario2.find(*it);
      if (it2==diccionario2.end()){
        it2 = diccionario.find(*it);
      }

      pareja sum;
			sum.first=0;
      int mayor=0;
      for (int k=0; k<V.SizeDomain(V.Consecuente()); k++){
        if (it2->second[k]>it2->second[mayor])
          mayor=k;
      }

      for (int k=0; k<V.SizeDomain(V.Consecuente()); k++){
        sum.first += it2->second[k];
      }
			sum.second=it2->second[mayor];
			sum.first = sum.first / sum.second;

      reves.insert(pair<pareja, unordered_map<string,vector<int> >::iterator >  (sum, it2) );
    }

    rFinal.clear();
    map<string, vector<double> > rFinalConPeso;

    int t = 1 ;
    vector<bool> cubierto;

    for (int i=0; i<Es.N_Examples(); i++){
      cubierto.push_back(false);
    }

    double porcentaje2;
    for (auto it=reves.rbegin(); it!=reves.rend(); it++){
      porcentaje2 = 100.0*t / reves.size();

      cout << porcentaje2 << " %\t"
           << InterpretaRegla(V,it->second->first) << "\t"
           << "\t" <<it->second->second[0]
           << "\t" << it->second->second[1]
           << endl;

      if (CubreEjemplos(it->second->first,V,Es,Tabla, cubierto)){
				double sum=0;
				int mayor=0;
	      for (int k=0; k<V.SizeDomain(V.Consecuente()); k++){
					sum += it->second->second[k];
	        if (it->second->second[k]>it->second->second[mayor])
	          mayor=k;
	      }

				vector<double> v;
				//Meto la clase
				v.push_back(mayor);

				//Meto ejemplos positivos
				v.push_back(it->second->second[mayor]);

				//Meto ejemplos negativos
				sum = sum - it->second->second[mayor];
				v.push_back(sum);
        rFinalConPeso.insert(pair<string, vector<double> > (it->second->first, v) );

				rFinal.insert(it->second->first);
      }
      else{
        deleteLine();
      }
      t++;
    }


    for (auto it=rFinal.begin(); it!=rFinal.end(); it++){
      cout << InterpretaRegla(V,*it) << endl;
    }

    return rFinalConPeso;

  }

  //-----------------------------------------------------------------------------------------------------

void GeneralizarPatrones(const example_set &Es, const VectorVar &V, double ***Tabla,
                         const multimap<pareja, unordered_map<string,vector<int> >::iterator, ComparaMultiMap> & reves,
                         const unordered_map<string,vector<int> > & diccionario,
                         unordered_map<string,vector<int> > & diccionario2,
                         set<string> &rFinal){
  double cero, uno, ncero, nuno;
	clock_t tiempo0, tiempo1;
	double tconsumido;
  string tmejor, smejor;

  double *varlevel = new double[V.N_Antecedente()+1];

  double porcentaje;
  int previoporcentaje;

  auto it2 = reves.rbegin();
  cout << endl;
  for (int t=0; t<reves.size(); t++){
    ncero =0;
    nuno = 0;
    //if (it2->second->second[1] > it2->second->second[0] ){
    if (true){

    porcentaje = 100.0* t / reves.size();
    previoporcentaje = porcentaje*100;
    porcentaje = previoporcentaje/100.0;

    string regla = it2->second->first;
    deleteLine();
    cout << porcentaje << " %\t"
         << InterpretaRegla(V,regla) << "\t" << endl;

    if (!IncluidoOnSet(rFinal, it2->second->first)){
			tconsumido=0;
      //it2->second->second[V.Consecuente()]=1;
      smejor = regla;


      for (int j=0; j<V.N_Antecedente(); j++)
        varlevel[j]=1;
      varlevel[V.N_Antecedente()]=0.5;

      int better = it2->second->second[0];
      double p_acierto = it2->second->second[0] / (1.0* it2->second->second[0] + it2->second->second[1]);
      double p_aciertoAux = it2->second->second[1] / (1.0* it2->second->second[0] + it2->second->second[1]);

      if (p_aciertoAux>p_acierto){
        p_acierto = p_aciertoAux;
        better = it2->second->second[1];
      }

      int  mejor;

      bool cambio=true, section1, section2;
      vector<int> tama;
      int longitud=0;
      for (int j=0; j<V.N_Antecedente(); j++){
        tama.push_back(longitud);
        longitud+=V.SizeDomain(j);
      }


      while (cambio){
        string regla2 = smejor;
        cambio = false;

        vector<int> orden;
        for (int j=0; j<V.N_Antecedente(); j++){
          orden.push_back(j);
        }
        std::random_shuffle(orden.begin(),orden.end());

        //#pragma omp parallel for default (none)
        for (int t=0; t<V.N_Antecedente() /*and !cambio*/; t++){
          int j=orden[t];
          //int j=t;
          if (varlevel[j]==1){
            int clase0=ncero, clase1=nuno;

            for (int l=0; l<V.SizeDomain(j); l++)
              regla2[tama[j]+l]='1';

            //cout << "->" << regla2 << endl;
            auto it4 = diccionario2.find(regla2);
            if (it4==diccionario2.end()){
              tconsumido += CalculoReglas(Es,V, regla2, Tabla, j, varlevel, diccionario, diccionario2, clase0, clase1);
            }
            else {
              //cout << "\t\t\tYa estaba\n ";
              clase0 = it4->second[0];
              clase1 = it4->second[1];
            }

            double p_aciertoNew = clase0 / (1.0*clase0+clase1);
            double p_aciertoNewAux = clase1 / (1.0*clase0+clase1);
            int magnitud = clase0;

            if (p_aciertoNewAux> p_aciertoNew){
              p_aciertoNew = p_aciertoNewAux;
              magnitud = clase1;
            }

            if (p_aciertoNew>p_acierto or (p_aciertoNew == p_acierto and magnitud>=better)){
              mejor = j;
              better = magnitud;
              p_acierto = p_aciertoNew;
              cambio=true;
              cero = clase0;
              uno = clase1;
              tmejor=regla2;
            }
          }
          regla2=smejor;

        }
        if (cambio){
          //cambio=false;
          varlevel[mejor]=0;
          smejor=tmejor;
          ncero=cero;
          nuno=uno;
          string aux;
          if (section1)
            aux = "Conjuntos ";
          else
            aux = "Inferencia";
          //cout << endl << "El mejor es j= " << mejor << endl;
          //Pausa();
          deleteLine();
          cout << porcentaje << " %\t"
               << InterpretaRegla(V,regla) << "\t"
               << "\t" << InterpretaRegla(V,smejor)
               << "\t" << cero
               << "\t" << uno
               << "\t" << tconsumido
               << "    \t" << "--"
               << "\t" << "--"
               << "\t" << "--"
               << "\t" << diccionario2.size()
               << "\t" << aux
               << endl;

        }
      }

      //Pausa();
      string regla2=regla;
      for (int j=0; j<V.N_Antecedente(); j++){
        if (varlevel[j]==0)
          for (int l=0; l<V.SizeDomain(j); l++)
            regla2[tama[j]+l]='1';
      }

      // Inserto la nueva regla
      rFinal.insert(regla2);

      if (rFinal.size()%10==0){
        SalvarMap(diccionario2);
        SalvarSet(rFinal);
      }

      //int inCero, inUno;
      //bool salir=false;
      //Infer(Es, V, Tabla, regla2, varlevel, inCero, inUno, salir);
      deleteLine();
      cout << porcentaje << " %\t"
           << InterpretaRegla(V,regla) << "\t"
           << "\t" << InterpretaRegla(V,smejor)
           << "\t" << cero
           << "\t" << uno
           << "\t" << tconsumido
           << "    \t" << rFinal.size()
           << "\t" << "--"
           << "\t" << "--"
           //<< "\t" << inCero
           //<< "\t" << inUno
           << "\t" << diccionario2.size()
           << endl;
    }
    else {
      deleteLine();
      //cout << InterpretaRegla(V,regla) << "\t\t" << "Ya Computado" << endl;
      //cout << porcentaje << " %\t"
      //     << InterpretaRegla(V,regla) << "\t" << endl;
    }
    //Pausa();
    cout << endl;
  }
  it2++;
  }

}

//-----------------------------------------------------------------------------------------------------

void PasarPatronesAReglas(const VectorVar &V, genetcode x, const map<string, vector<double> > & rFinalCP, ruleset &R){

	// Fijo el nivel de variable
  for (int i=0; i<V.N_Antecedente(); i++){
		x.PutValueReal(0,i,1);
	}
	x.PutValueReal(0,V.N_Antecedente(),0.5);

	// Nivel de Funcion
	int lon = x.TamaInteger(1);
	for (int i=0; i<lon; i++)
	  x.PutValueInteger(1,i,0);

	int tama = rFinalCP.begin()->first.size();

  for (auto it=rFinalCP.begin(); it!=rFinalCP.end(); it++){
		// Pongo la regla
		for (int i=0; i<tama; i++){
			x.PutValueBinary(0,i,it->first[i]);
		}
		// Pongo la clase
		int clase = it->second[0];
		x.PutValueInteger(0,0,clase);

		R.Add(x, it->second[1], it->second[2], false);
	}
}

//-----------------------------------------------------------------------------------------------------

set<string> PasarBasicosArFinal(const unordered_map<string,vector<int> > & diccionario){
	set <string> rFinal;
	for (auto it=diccionario.begin(); it!=diccionario.end(); it++){
		rFinal.insert(it->first);
	}
	return rFinal;
}


//-----------------------------------------------------------------------------------------------------


map<string, vector <double> >  ProcesarPatrones(const example_set &E, const example_set &ETest, const VectorVar &V, double *Contador, unordered_map <string, vector<int> > & diccionario){
	//example_set Es;
	//unordered_map<string,vector<int> > diccionario;
	//diccionario.max_load_factor(3.0);


 for (int i=0; i<15; i++)
   Contador[i]=0;


	//   for (int par=0; par<10; par++){

  Contador[11] = E.N_Examples();
	//Es = E.Extract_Training_Set(0,V.Consecuente(),0);
	//cout << "Creando la tabla de adaptaciones.....\n";
	double ***Tabla = Crear_Tabla_Adaptaciones(V,E);
	//double ***Tabla;
	//cout << "Terminada la tabla de adaptaciones.....\n";

	SacarPatronesBasicos(E, V, Tabla, diccionario);
	Destruir_Tabla_Adaptaciones(V,E, Tabla);
	//cout << "Numero de patrones: " << diccionario.size() << endl;
	Contador[0]=diccionario.size();

	multimap <pareja, unordered_map<string,vector<int> >::iterator, ComparaMultiMap> reves;
	OrdenarPatrones(V, diccionario, reves);

	auto it=reves.begin();
	int max=0;
	for (it=reves.begin(); it!=reves.end(); it++){
		double aux1=0, aux2=0, best=0;
		for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
			aux1=aux1+ it->second->second[i];
			if (it->second->second[i]>0){
				aux2++;
			}
			if (it->second->second[i]>it->second->second[best]){
				best=i;
			}
		}

	  if (aux1 == 1){
			Contador[3]++;
		}
		if (aux1> max){
			max = aux1;
		}
		if (aux2==1){
			Contador[1]++;
			Contador[4]+=aux1;
		}

		if (it->second->second[0] > it->second->second[1]){
			Contador[5]+=it->second->second[best];
			Contador[6]++;
		}
		else if (it->second->second[1] > it->second->second[0]){
			Contador[5]+=it->second->second[best];
			Contador[7]++;

		}
		else{
			Contador[5]+=it->second->second[best];
			Contador[8]++;
		}


		Contador[2]+=aux1;
	}

	if (Contador[0] == 0){
	  Contador[6]= 0;
		Contador[7]= 0;
		Contador[8]= 0;
	}
	else {
	  Contador[6]= 100.0*Contador[6]/Contador[0];
		Contador[7]= 100.0*Contador[7]/Contador[0];
		Contador[8]= 100.0*Contador[8]/Contador[0];
	}

  //cout << "Patrones con dominio de la clase 0: " << Contador[6] << endl;
	//cout << "Patrones con empate entre clases  : " << Contador[8] << endl;
	//cout << "Patrones con dominio de la clase 1: " << Contador[7] << endl;

	it = reves.end();
	it--;
	int aux1=0;
	for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
		aux1=aux1+ it->second->second[i];
	}
	Contador[9]= aux1;
	//cout << "Max.Red: " << Contador[9] << endl;

  it = reves.begin();
	for (int i=0; i<diccionario.size()/2; i++){
		it++;
	}

	aux1=0;
	for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
		aux1=aux1+ it->second->second[i];
	}
  Contador[10]= aux1;
	//cout << "Mediana Red: " << Contador[10] << endl;

	//cout << "Media por patron: " << Contador[2]/diccionario.size() << endl;

  Contador[3]= (100.0*Contador[3])/diccionario.size();
	Contador[1]= (100.0*Contador[1])/diccionario.size();
	Contador[4]= (100.0*Contador[4]);
	Contador[5]= (100.0*Contador[5]);
	//cout << "Pat_1_ejemplo: " << Contador[3] << endl;
	//cout << "Patrones 1Sola Clase: " << Contador[1] << endl;
	//cout << "Aciertos con clases exclusivas: " << Contador[4] << endl;
	//cout << "Aciertos tomando clase Mayoritaria: " << Contador[5] << endl;

	double ***TablaTes = Crear_Tabla_Adaptaciones(V,ETest);
  int no_cubiertos;
	Contador[14] = SacarPatronesBasicosTest(ETest, V, TablaTes, diccionario, no_cubiertos);
	Contador[13] = 100.0*Contador[14]/ETest.N_Examples();
	Contador[12] = 100.0*no_cubiertos/ETest.N_Examples();
  if (ETest.N_Examples()-no_cubiertos == 0){
		Contador[14] = 100;
	}
	else {
		Contador[14] = 100.0 *Contador[14] / (ETest.N_Examples() - no_cubiertos);
	}
	//cout << "Acierto sobre test: " << Contador[13] << endl;
	//cout << "Porcent nuevos patrones en test: " << Contador[12] << endl;
	Destruir_Tabla_Adaptaciones(V,ETest, TablaTes);

	//cout << "Terminado 1\n";

  //Pausa();

	/*set<string> rFinal;

	if (Existe("SetFile.dat")){
		diccionario2 = CargarMap();
		rFinal = CargarSet();
		//ReducirSet (Es,V, Tabla, diccionario, diccionario2, rFinal);
	}

	GeneralizarPatrones(Es, V, Tabla, reves, diccionario, diccionario2, rFinal);
	//rFinal = PasarBasicosArFinal(diccionario);

	cout << "Numero de reglas: " << rFinal.size() << endl;
	SalvarMap(diccionario2);
	SalvarSet(rFinal);

	// Reducir
	map <string, vector<double> > rFinalConPeso;
	rFinalConPeso = ReducirSet (Es,V, Tabla, diccionario, diccionario2, rFinal);
	SalvarSet(rFinal);

	Destruir_Tabla_Adaptaciones(V,Es, Tabla);
	Destruir_Tabla_Adaptaciones(V,ETest, TablaTes);*/

	map <string, vector<double> > rFinalConPeso;
	return rFinalConPeso;
}



map<string, vector <double> >  ProcesarPatronesTraining(const example_set &E,  const VectorVar &V, double *Contador, unordered_map <string, vector<int> > & diccionario){
	//example_set Es;
	//unordered_map<string,vector<int> > diccionario;
	//diccionario.max_load_factor(3.0);


 for (int i=0; i<15; i++)
   Contador[i]=0;


	//   for (int par=0; par<10; par++){

  Contador[11] = E.N_Examples();
	//Es = E.Extract_Training_Set(0,V.Consecuente(),0);
	//cout << "Creando la tabla de adaptaciones.....\n";
	double ***Tabla = Crear_Tabla_Adaptaciones(V,E);
	//double ***Tabla;
	//cout << "Terminada la tabla de adaptaciones.....\n";

	// Antes de aprender los nuevos patrones
	int no_cubiertos;
	Contador[14] = SacarPatronesBasicosTest(E, V, Tabla, diccionario, no_cubiertos);
	Contador[13] = 100.0*Contador[14]/E.N_Examples();
	Contador[12] = 100.0*no_cubiertos/E.N_Examples();
  if (E.N_Examples()-no_cubiertos == 0){
		Contador[14] = 100;
	}
	else {
		Contador[14] = 100.0 *Contador[14] / (E.N_Examples() - no_cubiertos);
	}


	// Aprender los nuevos patrones
	SacarPatronesBasicos(E, V, Tabla, diccionario);
	Destruir_Tabla_Adaptaciones(V,E, Tabla);
	//cout << "Numero de patrones: " << diccionario.size() << endl;
	Contador[0]=diccionario.size();

	multimap <pareja, unordered_map<string,vector<int> >::iterator, ComparaMultiMap> reves;
	OrdenarPatrones(V, diccionario, reves);

	auto it=reves.begin();
	int max=0;
	for (it=reves.begin(); it!=reves.end(); it++){
		double aux1=0, aux2=0, best=0;
		for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
			aux1=aux1+ it->second->second[i];
			if (it->second->second[i]>0){
				aux2++;
			}
			if (it->second->second[i]>it->second->second[best]){
				best=i;
			}
		}

	  if (aux1 == 1){
			Contador[3]++;
		}
		if (aux1> max){
			max = aux1;
		}
		if (aux2==1){
			Contador[1]++;
			Contador[4]+=aux1;
		}

		if (it->second->second[0] > it->second->second[1]){
			Contador[5]+=it->second->second[best];
			Contador[6]++;
		}
		else if (it->second->second[1] > it->second->second[0]){
			Contador[5]+=it->second->second[best];
			Contador[7]++;

		}
		else{
			Contador[5]+=it->second->second[best];
			Contador[8]++;
		}


		Contador[2]+=aux1;
	}

	if (Contador[0] == 0){
	  Contador[6]= 0;
		Contador[7]= 0;
		Contador[8]= 0;
	}
	else {
	  Contador[6]= 100.0*Contador[6]/Contador[0];
		Contador[7]= 100.0*Contador[7]/Contador[0];
		Contador[8]= 100.0*Contador[8]/Contador[0];
	}

  //cout << "Patrones con dominio de la clase 0: " << Contador[6] << endl;
	//cout << "Patrones con empate entre clases  : " << Contador[8] << endl;
	//cout << "Patrones con dominio de la clase 1: " << Contador[7] << endl;

	it = reves.end();
	it--;
	int aux1=0;
	for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
		aux1=aux1+ it->second->second[i];
	}
	Contador[9]= aux1;
	//cout << "Max.Red: " << Contador[9] << endl;

  it = reves.begin();
	for (int i=0; i<diccionario.size()/2; i++){
		it++;
	}

	aux1=0;
	for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
		aux1=aux1+ it->second->second[i];
	}
  Contador[10]= aux1;
	//cout << "Mediana Red: " << Contador[10] << endl;

	//cout << "Media por patron: " << Contador[2]/diccionario.size() << endl;

  Contador[3]= (100.0*Contador[3])/diccionario.size();
	Contador[1]= (100.0*Contador[1])/diccionario.size();
	Contador[4]= (100.0*Contador[4]);
	Contador[5]= (100.0*Contador[5]);
	//cout << "Pat_1_ejemplo: " << Contador[3] << endl;
	//cout << "Patrones 1Sola Clase: " << Contador[1] << endl;
	//cout << "Aciertos con clases exclusivas: " << Contador[4] << endl;
	//cout << "Aciertos tomando clase Mayoritaria: " << Contador[5] << endl;

	//cout << "Acierto sobre test: " << Contador[13] << endl;
	//cout << "Porcent nuevos patrones en test: " << Contador[12] << endl;

	//cout << "Terminado 1\n";

  //Pausa();

	/*set<string> rFinal;

	if (Existe("SetFile.dat")){
		diccionario2 = CargarMap();
		rFinal = CargarSet();
		//ReducirSet (Es,V, Tabla, diccionario, diccionario2, rFinal);
	}

	GeneralizarPatrones(Es, V, Tabla, reves, diccionario, diccionario2, rFinal);
	//rFinal = PasarBasicosArFinal(diccionario);

	cout << "Numero de reglas: " << rFinal.size() << endl;
	SalvarMap(diccionario2);
	SalvarSet(rFinal);

	// Reducir
	map <string, vector<double> > rFinalConPeso;
	rFinalConPeso = ReducirSet (Es,V, Tabla, diccionario, diccionario2, rFinal);
	SalvarSet(rFinal);

	Destruir_Tabla_Adaptaciones(V,Es, Tabla);
	Destruir_Tabla_Adaptaciones(V,ETest, TablaTes);*/

	map <string, vector<double> > rFinalConPeso;
	return rFinalConPeso;
}




void  ProcesarPatronesTest(const example_set &E,  const VectorVar &V, double *Contador, unordered_map <string, vector<int> > & diccionario){
	//example_set Es;
	//unordered_map<string,vector<int> > diccionario;
	//diccionario.max_load_factor(3.0);


 for (int i=12; i<15; i++)
   Contador[i]=0;



	double ***TablaTes = Crear_Tabla_Adaptaciones(V,E);
  int no_cubiertos;
	Contador[14] = SacarPatronesBasicosTest(E, V, TablaTes, diccionario, no_cubiertos);
	Contador[13] = 100.0*Contador[14]/E.N_Examples();
	Contador[12] = 100.0*no_cubiertos/E.N_Examples();
  if (E.N_Examples()-no_cubiertos == 0){
		Contador[14] = 100;
	}
	else {
		Contador[14] = 100.0 *Contador[14] / (E.N_Examples() - no_cubiertos);
	}
	//cout << "Acierto sobre test: " << Contador[13] << endl;
	//cout << "Porcent nuevos patrones en test: " << Contador[12] << endl;
	Destruir_Tabla_Adaptaciones(V,E, TablaTes);

}


//-----------------------------------------------------------------------------------------------------
double Average (double *v, int n_datos){
	double suma = 0;
	for (int i=0; i<n_datos; i++){
		suma = suma + v[i];
	}

	return suma / n_datos;
}


//-----------------------------------------------------------------------------------------------------
double SDesviation (double *v, int n_datos){
	double suma = 0;
	double media = Average(v,n_datos);

	for (int i=0; i<n_datos; i++){
		suma = suma + pow(v[i]-media, 2);
	}

	return sqrt(suma / n_datos);
}

//-----------------------------------------------------------------------------------------------------

double Average (int *v, int n_datos){
	double suma = 0;
	for (int i=0; i<n_datos; i++){
		suma = suma + v[i];
	}

	return suma / n_datos;
}


//-----------------------------------------------------------------------------------------------------
double SDesviation (int *v, int n_datos){
	double suma = 0;
	double media = Average(v,n_datos);

	for (int i=0; i<n_datos; i++){
		suma = suma + pow(v[i]-media, 2);
	}

	return sqrt(suma / n_datos);
}


//-----------------------------------------------------------------------------------------------------

void CalculoEstimadoSobrePatrones(example_set & E, const VectorVar &V, int num_par){
	unordered_map<string,vector<int> > diccionario;
	diccionario.max_load_factor(3.0);
	double reSult[15], SumReSult[15]={0};
	for (int par=0; par<num_par; par++){
		cout <<"\n===============\n Ejecucion " << par << "\n===============\n";

		example_set E_Par_Completo=E.Extract_Training_Set(par, V.Consecuente(), 0);
		example_set E_Par_Test = E.Extract_Test_Set(par, V.Consecuente());
		for (int i=0; i<15; i++){
			reSult[i]=0;
		}
		// Parte para ver como evolucionan la recepcion de patrones
		ProcesarPatronesTraining(E_Par_Completo, V, reSult, diccionario);
		ProcesarPatronesTest(E_Par_Test, V, reSult, diccionario);
		reSult[11] = E_Par_Completo.N_Examples();
		cerr << "Numero medio de Ejemplos: " << reSult[11] << endl;
		cerr << "Numero medio de Patrones: " << reSult[0] << endl;
		cerr << "         Mayoria Clase 0: " << reSult[6] << endl;
		cerr << "                  Empate: " << reSult[8] << endl;
		cerr << "         Mayoria Clase 1: " << reSult[7] << endl;
		cerr << "Maxima Redundancia      : " << reSult[9] << endl;
		cerr << "Media ejemplos/patron   : " << reSult[11]/reSult[0] << endl;
		cerr << "Mediana Redundancia     : " << reSult[10] << endl;
		cerr << "Patrones con 1 ejemplo  : " << reSult[3] << endl;
		cerr << "Patrones de 1 sola clase: " << reSult[1] << endl;
		cerr << "Media Train sólo 1 clase: " << reSult[4] << endl;
		cerr << "Media Train clase mayor : " << reSult[5] << endl;
		cerr << "% Acierto Test (Global) : " << reSult[13] << endl;
		cerr << "% nuevos patrones Test  : " << reSult[12] << endl;
		cerr << "% Test sin nocubiertos  : " << reSult[14] << endl;
		cerr << "%Error intrinseco       : " << 100 -reSult[14] << endl;
		cerr << endl;
		//Pausa();


		for (int i=0; i<15; i++){
			SumReSult[i]+=reSult[i];
			reSult[i] = 0;
		}
	}

	for (int i=0; i<15; i++){
		SumReSult[i]= 1.0*SumReSult[i]/10;
	}
	cerr << endl;

	cerr << "Numero medio de Ejemplos: " << SumReSult[11] << endl;
	cerr << "Numero medio de Patrones: " << SumReSult[0] << endl;
	cerr << "         Mayoria Clase 0: " << SumReSult[6] << endl;
	cerr << "                  Empate: " << SumReSult[8] << endl;
	cerr << "         Mayoria Clase 1: " << SumReSult[7] << endl;
	cerr << "Maxima Redundancia      : " << SumReSult[9] << endl;
	cerr << "Media ejemplos/patron   : " << SumReSult[11]/SumReSult[0] << endl;
	cerr << "Mediana Redundancia     : " << SumReSult[10] << endl;
	cerr << "Patrones con 1 ejemplo  : " << SumReSult[3] << endl;
	cerr << "Patrones de 1 sola clase: " << SumReSult[1] << endl;
	cerr << "Media Train sólo 1 clase: " << SumReSult[4] << endl;
	cerr << "Media Train clase mayor : " << SumReSult[5] << endl;
	cerr << "% Acierto Test (Global) : " << SumReSult[13] << endl;
	cerr << "% nuevos patrones Test  : " << SumReSult[12] << endl;
	cerr << "% Test sin nocubiertos  : " << SumReSult[14] << endl;
	cerr << "%Error intrinseco       : " << 100 -SumReSult[14] << endl;
	cerr << endl;



	//Pausa();

	//exit(0);

}
//-----------------------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------------------

vector<vector<double> > CalculoEstimadoSobrePatrones2(example_set & E, const VectorVar &V, int num_par){
  vector<TestResult> result;
	TestResult porDefecto, aux;
	Pattern P(V);

	vector <vector <double> > salida;
	vector<double> s;

  s.clear();
	for (int i=0; i<V.SizeDomain(V.Consecuente())+1; i++){
		salida.push_back(s);
	}

	/*P.ExtraerPatronesBasicos(E,V, porDefecto);
	cerr << "Numero de Ejemplos: " << E.N_Examples() << endl;
	for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
		cerr << "           clase " << i << ": " << P.N_Ejemplos_paraClase(i) << " (" <<100.0*P.N_Ejemplos_paraClase(i) / E.N_Examples() << ")" <<endl;
		salida[i+1].push_back(P.N_Ejemplos_paraClase(i));
	}
	cerr << "Numero de Patrones: " << P.N_Pattern() << endl;
	for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
		cerr << "           clase " << i << ": " << P.N_Patrones_paraClase(i) << " (" <<100.0*P.N_Patrones_paraClase(i) / P.N_Pattern() << ")" <<endl;
		salida[i+1].push_back(P.N_Patrones_paraClase(i));
	}
	cerr << "Media ejemplos/patron   : " << 1.0*E.N_Examples()/P.N_Pattern() << endl;
	//P.TestearPatronesBasicos(E,V,porDefecto);
	PintaResultadosTraining(porDefecto);

	for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
		salida[i+1].push_back(100.0 * porDefecto.acc[i]/(porDefecto.cubiertos[i]));
		salida[i+1].push_back(100.0 - 100.0 * porDefecto.acc[i]/(porDefecto.cubiertos[i]));
	}
  */
	//Pausa();

  clock_t crono, etapa;


  crono = 0;
	for (int par=0; par<num_par; par++){

		Pattern Patrones(V);
		cout <<"\n===============\n Ejecucion " << par << "\n===============\n";

		example_set E_Par_Completo=E.Extract_Training_Set(par, V.Consecuente(), 0);
		example_set E_Par_Test = E.Extract_Test_Set(par, V.Consecuente());

		result.push_back(porDefecto);

		etapa = clock();
		Patrones.ExtraerPatronesBasicos(E_Par_Completo,V, result[par]);
		etapa = clock() - etapa;
		crono += etapa;
		//Patrones.TestearPatronesBasicos(E_Par_Completo,V,result[par]);

		cerr << "Numero de Ejemplos: " << E_Par_Completo.N_Examples() << endl;
		for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
			cerr << "         clase " << i << ": " << Patrones.N_Ejemplos_paraClase(i) << " (" <<100.0*Patrones.N_Ejemplos_paraClase(i) / E_Par_Completo.N_Examples() << ")" <<endl;
		}
		cerr << "Numero de Patrones: " << Patrones.N_Pattern() << endl;
		for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
			cerr << "         clase " << i << ": " << Patrones.N_Patrones_paraClase(i) << " (" <<100.0*Patrones.N_Patrones_paraClase(i) / Patrones.N_Pattern() << ")" <<endl;
		}
		cerr << "Media ejemplos/patron   : " << 1.0*E_Par_Completo.N_Examples()/Patrones.N_Pattern() << endl;
		for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
			cerr << "         clase " << i << ": " << 1.0 * Patrones.N_Ejemplos_paraClase(i) / Patrones.N_Patrones_paraClase(i) <<endl;
		}
		cerr << "Porcentaje sobre el global de patrones: " << 100.0 * Patrones.N_Pattern() / P.N_Pattern() << endl;
		for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
			cerr << "         clase " << i << ": " << Patrones.N_Patrones_paraClase(i) << " (" <<100.0*Patrones.N_Patrones_paraClase(i) / P.N_Pattern() << ")" <<endl;
		}
		//cerr << "Media Train clase mayor : " << result[par].acierto_global << endl;
		cerr << "\t-------------------- Training -----------" << endl;
		PintaResultadosTraining(result[par]);
		cerr << "Tiempo: " << 1.0*(etapa)/CLOCKS_PER_SEC << endl;

		cerr << "\n\t--------------------- Test  -----------" << endl;
		Patrones.TestearPatronesBasicos(E_Par_Test,V,result[par]);
		PintaResultadosTest(result[par], true);
		/*cerr << "% Acierto Test (Global) : " << result[par].acierto_global << endl;
		cerr << "% nuevos patrones Test  : " << result[par].porcentaje_nuevos_patrones << endl;
		cerr << "% Test sin nocubiertos  : " << result[par].acierto_sinNoCubiertos << endl;
		cerr << "%Error intrinseco       : " << result[par].error_intrinseco << endl;*/
		cerr << endl;
		//Pausa();
}

	/*cerr << endl;

	cerr << "Numero medio de Ejemplos: " << SumReSult[11] << endl;
	cerr << "Numero medio de Patrones: " << SumReSult[0] << endl;
	cerr << "         Mayoria Clase 0: " << SumReSult[6] << endl;
	cerr << "                  Empate: " << SumReSult[8] << endl;
	cerr << "         Mayoria Clase 1: " << SumReSult[7] << endl;
	cerr << "Maxima Redundancia      : " << SumReSult[9] << endl;
	cerr << "Media ejemplos/patron   : " << SumReSult[11]/SumReSult[0] << endl;
	cerr << "Mediana Redundancia     : " << SumReSult[10] << endl;
	cerr << "Patrones con 1 ejemplo  : " << SumReSult[3] << endl;
	cerr << "Patrones de 1 sola clase: " << SumReSult[1] << endl;
	cerr << "Media Train sólo 1 clase: " << SumReSult[4] << endl;
	cerr << "Media Train clase mayor : " << SumReSult[5] << endl;
	cerr << "% Acierto Test (Global) : " << SumReSult[13] << endl;
	cerr << "% nuevos patrones Test  : " << SumReSult[12] << endl;
	cerr << "% Test sin nocubiertos  : " << SumReSult[14] << endl;
	cerr << "%Error intrinseco       : " << 100 -SumReSult[14] << endl;
	cerr << endl;
*/

cerr << "\n\t------------------ Media Test -----------" << endl;
TestResult aux2;
cerr << "Tiempo Medio de Aprendizaje: " << (1.0*(crono)/CLOCKS_PER_SEC)/10 << endl;
CalcularMediaTestResult(result, aux2);
PintaResultadosTest(aux2, true);

s.clear();
s.push_back(E.N_Examples());
s.push_back(P.N_Pattern());
s.push_back(porDefecto.acierto_global);
s.push_back(porDefecto.error_intrinseco);
s.push_back(aux2.acierto_global);
s.push_back(aux2.acierto_sinNoCubiertos);
s.push_back(aux2.error_intrinseco);
s.push_back(aux2.porcentaje_nuevos_patrones);
salida[0] = s;

for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
	salida[i+1].push_back(100.0 * aux2.acc[i]/(aux2.cubiertos[i]+aux2.no_cubiertos[i]));
	if (aux2.cubiertos[i]>0){
		salida[i+1].push_back(100.0 * aux2.acc[i]/(aux2.cubiertos[i]));
		salida[i+1].push_back(100.0 - 100.0 * aux2.acc[i]/(aux2.cubiertos[i]));
	}
	else {
		salida[i+1].push_back(0);
		salida[i+1].push_back(100);
	}
	salida[i+1].push_back(100.0 * aux2.no_cubiertos[i]/(aux2.cubiertos[i]+aux2.no_cubiertos[i]));
}


//Pausa();

//exit(0);

return salida;

}

string PonerTiempo (double x){
	int t = x;
	int horas = t / 3600;
	t = t - horas * 3600;
	int minutos = t / 60;
	t = t - minutos * 60;
	int seg = t;
	string sseg, smin, shoras;
	if (seg<10){
	  sseg = "0"+to_string(seg);
  }
	else {
		sseg = to_string(seg);
	}

	if (minutos<10){
		smin = "0"+to_string(minutos);
	}
	else{
		smin = to_string(minutos);
	}


	string salida = to_string(horas) + ":" + smin + ":" + sseg;
	return salida;
}



vector<vector<double> > CalculoEstimadoSobrePatrones2_Para5subconjuntos(example_set & E, const VectorVar &V5, const VectorVar &V3, const VectorVar &V2, int num_par){
  vector<TestResult> result, result2, resultTrain, resultTrain2;
	TestResult porDefecto, aux;
	Pattern P(V5);

	vector <vector <double> > salida;
	vector<double> s;

  s.clear();
	for (int i=0; i<V5.SizeDomain(V5.Consecuente())+1; i++){
		salida.push_back(s);
	}

	/*P.ExtraerPatronesBasicos(E,V, porDefecto);
	cerr << "Numero de Ejemplos: " << E.N_Examples() << endl;
	for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
		cerr << "           clase " << i << ": " << P.N_Ejemplos_paraClase(i) << " (" <<100.0*P.N_Ejemplos_paraClase(i) / E.N_Examples() << ")" <<endl;
		salida[i+1].push_back(P.N_Ejemplos_paraClase(i));
	}
	cerr << "Numero de Patrones: " << P.N_Pattern() << endl;
	for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
		cerr << "           clase " << i << ": " << P.N_Patrones_paraClase(i) << " (" <<100.0*P.N_Patrones_paraClase(i) / P.N_Pattern() << ")" <<endl;
		salida[i+1].push_back(P.N_Patrones_paraClase(i));
	}
	cerr << "Media ejemplos/patron   : " << 1.0*E.N_Examples()/P.N_Pattern() << endl;
	//P.TestearPatronesBasicos(E,V,porDefecto);
	PintaResultadosTraining(porDefecto);

	for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
		salida[i+1].push_back(100.0 * porDefecto.acc[i]/(porDefecto.cubiertos[i]));
		salida[i+1].push_back(100.0 - 100.0 * porDefecto.acc[i]/(porDefecto.cubiertos[i]));
	}
  */
	//Pausa();

  clock_t crono, etapa;


  crono = 0;
	for (int par=0; par<num_par; par=par+2){

		Pattern Patrones5(V5), Patrones3(V3), Patrones2(V2);
		cout <<"\n===============\n Ejecucion " << par << "\n===============\n";

		example_set E_Par_Completo=E.Extract_Training_Set2(par,par+1, V5.Consecuente(), 0);
		example_set E_Par_Test = E.Extract_Test_Set2(par, par+1, V5.Consecuente());

		result.push_back(porDefecto);
		result2.push_back(porDefecto);
		resultTrain.push_back(porDefecto);
		resultTrain2.push_back(porDefecto);


		etapa = clock();
		cout << "...... Aprendiendo con 2 etiquetas .....\n";
		Patrones2.ExtraerPatronesBasicos(E_Par_Completo,V2, resultTrain[par/2]);
		cout << "...... Aprendiendo con 3 etiquetas .....\n";
		Patrones3.ExtraerPatronesBasicos(E_Par_Completo,V3, resultTrain[par/2]);
		cout << "...... Aprendiendo con 5 etiquetas .....\n";
		Patrones5.ExtraerPatronesBasicos(E_Par_Completo,V5, resultTrain[par/2]);
		etapa = clock() - etapa;
		crono += etapa;
		//Patrones.TestearPatronesBasicos(E_Par_Completo,V,result[par]);
		Patrones5.TestearRecursivoVariosDicionarios(E_Par_Completo,V5,Patrones3,V3,Patrones2, V2, resultTrain2[par/2]);

		cerr << "Numero de Ejemplos: " << E_Par_Completo.N_Examples() << endl;
		for (int i=0; i<V5.SizeDomain(V5.Consecuente()); i++){
			cerr << "         clase " << i << ": " << Patrones5.N_Ejemplos_paraClase(i) << " (" <<100.0*Patrones5.N_Ejemplos_paraClase(i) / E_Par_Completo.N_Examples() << ")" <<endl;
		}
		cerr << "Numero de Patrones: " << Patrones5.N_Pattern() << endl;
		for (int i=0; i<V5.SizeDomain(V5.Consecuente()); i++){
			cerr << "         clase " << i << ": " << Patrones5.N_Patrones_paraClase(i) << " (" <<100.0*Patrones5.N_Patrones_paraClase(i) / Patrones5.N_Pattern() << ")" <<endl;
		}
		cerr << "Media ejemplos/patron   : " << 1.0*E_Par_Completo.N_Examples()/Patrones5.N_Pattern() << endl;
		for (int i=0; i<V5.SizeDomain(V5.Consecuente()); i++){
			cerr << "         clase " << i << ": " << 1.0 * Patrones5.N_Ejemplos_paraClase(i) / Patrones5.N_Patrones_paraClase(i) <<endl;
		}
		cerr << "Porcentaje sobre el global de patrones: " << 100.0 * Patrones5.N_Pattern() / P.N_Pattern() << endl;
		for (int i=0; i<V5.SizeDomain(V5.Consecuente()); i++){
			cerr << "         clase " << i << ": " << Patrones5.N_Patrones_paraClase(i) << " (" <<100.0*Patrones5.N_Patrones_paraClase(i) / P.N_Pattern() << ")" <<endl;
		}
		//cerr << "Media Train clase mayor : " << result[par].acierto_global << endl;
		cerr << "\t-------------------- Training -----------" << endl;
		PintaResultadosTraining(resultTrain[par/2]);
		cerr << "Tiempo: " << 1.0*(etapa)/CLOCKS_PER_SEC << endl;

		cerr << "\n\t--------------------- Test  -----------" << endl;
		Patrones5.TestearPatronesBasicos(E_Par_Test,V5,result[par/2]);
		//Patrones2.TestearPatronesBasicos(E_Par_Test,V2,result[par/2]);
		PintaResultadosTest(result[par/2], true);
		Patrones5.TestearRecursivoVariosDicionarios(E_Par_Test,V5,Patrones3,V3,Patrones2,V2,result2[par/2]);
		PintaResultadosTest(result2[par/2], true);
		//Patrones.InferirDicionario(E_Par_Test,V);
		/*cerr << "% Acierto Test (Global) : " << result[par].acierto_global << endl;
		cerr << "% nuevos patrones Test  : " << result[par].porcentaje_nuevos_patrones << endl;
		cerr << "% Test sin nocubiertos  : " << result[par].acierto_sinNoCubiertos << endl;
		cerr << "%Error intrinseco       : " << result[par].error_intrinseco << endl;*/
		cerr << endl;
		//Pausa();
}

	/*cerr << endl;

	cerr << "Numero medio de Ejemplos: " << SumReSult[11] << endl;
	cerr << "Numero medio de Patrones: " << SumReSult[0] << endl;
	cerr << "         Mayoria Clase 0: " << SumReSult[6] << endl;
	cerr << "                  Empate: " << SumReSult[8] << endl;
	cerr << "         Mayoria Clase 1: " << SumReSult[7] << endl;
	cerr << "Maxima Redundancia      : " << SumReSult[9] << endl;
	cerr << "Media ejemplos/patron   : " << SumReSult[11]/SumReSult[0] << endl;
	cerr << "Mediana Redundancia     : " << SumReSult[10] << endl;
	cerr << "Patrones con 1 ejemplo  : " << SumReSult[3] << endl;
	cerr << "Patrones de 1 sola clase: " << SumReSult[1] << endl;
	cerr << "Media Train sólo 1 clase: " << SumReSult[4] << endl;
	cerr << "Media Train clase mayor : " << SumReSult[5] << endl;
	cerr << "% Acierto Test (Global) : " << SumReSult[13] << endl;
	cerr << "% nuevos patrones Test  : " << SumReSult[12] << endl;
	cerr << "% Test sin nocubiertos  : " << SumReSult[14] << endl;
	cerr << "%Error intrinseco       : " << 100 -SumReSult[14] << endl;
	cerr << endl;
*/

cerr << "\n\t------------------ Media Training 1 -----------" << endl;
TestResult aux2;
CalcularMediaTestResult(resultTrain, aux2);
PintaResultadosTest(aux2, true);


cerr << "\n\t------------------ Media Test 1 -----------" << endl;
cerr << "Tiempo Medio de Aprendizaje: " << (1.0*(crono)/CLOCKS_PER_SEC)/5 << endl;
CalcularMediaTestResult(result, aux2);
PintaResultadosTest(aux2, true);

cerr << "\n\t------------------ Media Training 2 -----------" << endl;
CalcularMediaTestResult(resultTrain2, aux2);
PintaResultadosTest(aux2, true);

cerr << "\n\t------------------ Media Test 2 -----------" << endl;
CalcularMediaTestResult(result2, aux2);
PintaResultadosTest(aux2, true);

s.clear();
s.push_back(E.N_Examples());
s.push_back(P.N_Pattern());
s.push_back(porDefecto.acierto_global);
s.push_back(porDefecto.error_intrinseco);
s.push_back(aux2.acierto_global);
s.push_back(aux2.acierto_sinNoCubiertos);
s.push_back(aux2.error_intrinseco);
s.push_back(aux2.porcentaje_nuevos_patrones);
salida[0] = s;

for (int i=0; i<V5.SizeDomain(V5.Consecuente()); i++){
	salida[i+1].push_back(100.0 * aux2.acc[i]/(aux2.cubiertos[i]+aux2.no_cubiertos[i]));
	if (aux2.cubiertos[i]>0){
		salida[i+1].push_back(100.0 * aux2.acc[i]/(aux2.cubiertos[i]));
		salida[i+1].push_back(100.0 - 100.0 * aux2.acc[i]/(aux2.cubiertos[i]));
	}
	else {
		salida[i+1].push_back(0);
		salida[i+1].push_back(100);
	}
	salida[i+1].push_back(100.0 * aux2.no_cubiertos[i]/(aux2.cubiertos[i]+aux2.no_cubiertos[i]));
}


//Pausa();

//exit(0);

return salida;

}



//-----------------------------------------------------------------------------------------------------

vector<vector<double> > CalculoEstimadoSobrePatronesConNEtiquetas(int nlabels, example_set & E, const VectorVar &V, int num_par){
  vector<TestResult> result, resultRecur;
	TestResult porDefecto, aux, porDefecto2;
	VectorVar V2;
	V2 = V;

	// Creo el dominio con n etiquetas en las variables continuas
	for (int i = 0; i<V.N_Antecedente(); i++){
		if (!V.Variable(i).IsDiscrete()){
			variable_t aux;
			aux.Asigna(nlabels, V.Variable(i).Inf_Range(), V.Variable(i).Sup_Range(), true, true, V.Variable(i).Name());
			V2.Asigna(i, aux);
		}
		else {
			V2.Asigna(i, V.Variable(i));
		}
	}

	int numeroTotalDePatrones;
	vector <vector <double> > salida;
	vector<double> s;

	s.clear();
	for (int i=0; i<V.SizeDomain(V.Consecuente())+1; i++){
		salida.push_back(s);
	}
	int redundancia = 1;

  {// Creo un ambito solo operar P2 y reducir el uso de la memoria automatica
		Pattern P2(V2);


		//P2.ExtraerPatronesBasicos(E,V2, porDefecto);
		//P2.ExtraerPatronesBasicosOriginalWM(E,V2, porDefecto);

		cout << "antes de ejecutar el aprendizaje TFM\n";
		P2.ExtraerPatronesBasicosAproximacionTFMRuben(E,V2, porDefecto, redundancia);

		numeroTotalDePatrones = P2.N_Pattern();
		cerr << "Numero de Ejemplos: " << E.N_Examples() << endl;
		for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
			cerr << "           clase " << i << ": " << P2.N_Ejemplos_paraClase(i) << " (" <<100.0*P2.N_Ejemplos_paraClase(i) / E.N_Examples() << ")" <<endl;
			salida[i+1].push_back(P2.N_Ejemplos_paraClase(i));
		}
		cerr << "Numero de Patrones: " << P2.N_Pattern() << endl;
		for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
			cerr << "           clase " << i << ": " << P2.N_Patrones_paraClase(i) << " (" <<100.0*P2.N_Patrones_paraClase(i) / P2.N_Pattern() << ")" <<endl;
			salida[i+1].push_back(P2.N_Patrones_paraClase(i));
		}
		cerr << "Media ejemplos/patron(2)   : " << 1.0*E.N_Examples()/P2.N_Pattern() << endl;
		PintaResultadosTraining(porDefecto);

		for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
			salida[i+1].push_back(100.0 * porDefecto.acc[i]/(porDefecto.cubiertos[i]));
			salida[i+1].push_back(100.0 - 100.0 * porDefecto.acc[i]/(porDefecto.cubiertos[i]));
		}
	}

	//Pausa();

  clock_t inicio, fin, timetrain = 0, timetestR = 0, timetestC = 0;
	double igual_clase=0, igual_adaptacion=0, igual_regla=0;
  double media_patrones = 0;


	for (int par=0; par<num_par; par++){
		Pattern Patrones2(V2);
		cout <<"\n===============\n Ejecucion " << par << "\n===============\n";

		example_set E_Par_Completo=E.Extract_Training_Set(par, V.Consecuente(), 0);
		example_set E_Par_Test = E.Extract_Test_Set(par, V.Consecuente());

		resultRecur.push_back(porDefecto);
		result.push_back(porDefecto);
		inicio = clock();
		//Patrones2.ExtraerPatronesBasicos(E_Par_Completo,V2, result[par]);
		//Patrones2.ExtraerPatronesBasicosOriginalWM(E_Par_Completo,V2, result[par]);

		cout << "antes de ejecutar el aprendizaje TFM\n";
		Patrones2.ExtraerPatronesBasicosAproximacionTFMRuben(E_Par_Completo,V2, resultRecur[par], redundancia);
		media_patrones += Patrones2.N_Pattern();


		fin = clock();
		timetrain += fin - inicio;

		cerr << "Numero de Ejemplos: " << E_Par_Completo.N_Examples() << endl;
		for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
			cerr << "         clase " << i << ": " << Patrones2.N_Ejemplos_paraClase(i) << " (" <<100.0*Patrones2.N_Ejemplos_paraClase(i) / E_Par_Completo.N_Examples() << ")" <<endl;
		}
		cerr << "Numero de Patrones: " << Patrones2.N_Pattern() << endl;
		for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
			cerr << "         clase " << i << ": " << Patrones2.N_Patrones_paraClase(i) << " (" <<100.0*Patrones2.N_Patrones_paraClase(i) / Patrones2.N_Pattern() << ")" <<endl;
		}
		cerr << "Media ejemplos/patron   : " << 1.0*E_Par_Completo.N_Examples()/(Patrones2.N_Pattern()) << endl;
		for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
			cerr << "         clase " << i << ": " << 1.0 * Patrones2.N_Ejemplos_paraClase(i) / Patrones2.N_Patrones_paraClase(i) <<endl;
		}
		cerr << "Porcentaje sobre el global de patrones: " << 100.0 * (Patrones2.N_Pattern()) / (numeroTotalDePatrones) << endl;
		for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
			cerr << "         clase " << i << ": " << Patrones2.N_Patrones_paraClase(i) << " (" <<100.0*Patrones2.N_Patrones_paraClase(i) / (numeroTotalDePatrones) << ")" <<endl;
		}
		//cerr << "Media Train clase mayor : " << result[par].acierto_global << endl;
		cerr << "\t-------------------- Training -----------" << endl;
		PintaResultadosTraining(result[par]);

		cerr << "\n\t--------------------- Test  -----------" << endl;
		double peso;
		resultRecur[par].cubiertos.clear();
	  resultRecur[par].no_cubiertos.clear();
	  resultRecur[par].acc.clear();
	  resultRecur[par].error_intrinseco_porClase.clear();

	  for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
	    resultRecur[par].cubiertos.push_back(0);
	    resultRecur[par].no_cubiertos.push_back(0);
	    resultRecur[par].acc.push_back(0);
	    resultRecur[par].error_intrinseco_porClase.push_back(0);
	  }

		// Calcular el peso y la clase de cada regla de las aprendidas para el inferencia CLASICA
		vector<double> pesoR;
		vector<int> claseR;
		Patrones2.CalcularPesoYClases(pesoR,claseR,true);

		clock_t t_salidas = 0;
		vector<clock_t> cronometro;

    inicio = clock();

    clock_t inicio_trozos = 0;
		clock_t step, step_begin = 0, step_end = 0;
		//Descomentar desde aqui hasta (A)
		vector<pair<int,pair <string, double> > >DisparosRecursiva;

		int partes = 1;
		double percentTestSet = 0.001;  //Porcentaje sobre el total del test sobre el que se está experimentando
		cout << "Se va a testear sobre el " << 100*percentTestSet*partes << " por ciento del total de ejemplos. Son " << (int) (E_Par_Test.N_Examples()*percentTestSet) << " ejemplos.\n\n";
		long int nex = E_Par_Test.N_Examples()*percentTestSet;
		long int old_porcion = 0;

		for (int t=0; t<partes; t++){ // Numero de trozos
			clock_t timePorcion = clock();
			for (int i=t*nex; i<(t+1)*nex; i++){
				// Opcion Testeo Original
				//int cl = Patrones2.TestearUnEjemploPatronesBasicos(E_Par_Test, i, V2, peso);

				// Opcion Testeo Recursivo
			  long int porcion = (10000.0*i/ ((t+1)*nex));
				double portion = porcion/100.0;
				clock_t projeccion=0;
				double t_uno = 0;
				if (i > 0){
					t_uno = 1.0* (inicio_trozos)/i;
					projeccion = t_uno * ((t+1)*nex-i);
				}

				if (old_porcion != porcion){
					clock_t b = clock();
					deleteLine();
					cout << "\tEx " << portion
					     << "\%  ETA: " << PonerTiempo(projeccion/CLOCKS_PER_SEC)
							 << "   Uno: " << t_uno / CLOCKS_PER_SEC
							 << "   TiempoTotal: " <<  PonerTiempo(t_uno * (t+1)*nex/CLOCKS_PER_SEC) << endl;
					old_porcion = porcion;
					t_salidas += (b - clock());
				}

				string regla="";
		    double umbral = 0.0;
				int NumMaxEvaluatedRules=1024;  // Numero de reglas de la base de reglas encontrada
		    int clase_predicha = -1;
				bool PCF = true;
				double grado = 0;
				string antecedenteSeleccionado = "-";
				step_begin = clock();


				//------------------------ Version Beam ------------
				int DistanciaHamming = 2;
				if (clase_predicha == -1)
				clase_predicha = Patrones2.InferenciaRecursivaOptimalizadaConProfundidadLimitada(E_Par_Test,V2,i,PCF,grado,antecedenteSeleccionado,DistanciaHamming);
				//------------------------ TFM Ruben ---------------
				if (clase_predicha == -1)
				clase_predicha = Patrones2.InferenciaRecursivaOptimalizada(E_Par_Test,V2,i,PCF,grado,antecedenteSeleccionado, NumMaxEvaluatedRules);

				//------------------------ Inferencia clasica ------
				if (clase_predicha == -1)
  			clase_predicha = Patrones2.TestearPatronesBasicosClassic_UnEjemplo(E_Par_Test,V2,i,pesoR,claseR, grado,antecedenteSeleccionado);

				//------------------------ Congreso ----------------
				//Patrones2.TestearRecursivoUnEjemplo(E_Par_Test, V2, i, regla, 0, 1, umbral, clase_predicha, PCF, antecedenteSeleccionado, NumMaxEvaluatedRules); grado = umbral;



        //------------------------ Prueba de sacar las salida por vecindarios
				//clase_predicha = Patrones2.InferenciaRecursivaOptimalizadaConProfundidadLimitada_SalidaPorProfundidades(E_Par_Test,V2,i,PCF,grado,antecedenteSeleccionado,DistanciaHamming);
				step_end = clock();
				step = step_end - step_begin;
				inicio_trozos += step;

				pair<int, pair<string, double> > item_disparo;
				int cl = clase_predicha;

				item_disparo.first = cl;
				item_disparo.second.first = antecedenteSeleccionado;
				item_disparo.second.second = grado;
				DisparosRecursiva.push_back(item_disparo);


				int clase = E_Par_Test.Data(i,V.Consecuente());
				if (cl == -1)	{
					resultRecur[par].no_cubiertos[clase]++;
				}
				else {
					resultRecur[par].cubiertos[clase]++;
					if (cl == clase){
						resultRecur[par].acc[clase]++;
					}
				}

			}  // Este es el final de (A) para la inferencia recursiva.
			timePorcion = 1.0*(clock() - timePorcion);
			ProcesarResultados(resultRecur[par]);
			fstream f;
			f.open("patrBasicos.csv", std::fstream::out | std::fstream::app);

			f << "\tPart: " << par << " Iter: " << t << " hasta: " << (t+1)*nex-1
		 	  << " Acc: "<< resultRecur[par].acierto_global << " AccSNC: " << resultRecur[par].acierto_sinNoCubiertos
			     << " err_int: " << resultRecur[par].error_intrinseco << " NewPattern: " << resultRecur[par].porcentaje_nuevos_patrones
					 << " tiempo: " << (1.0*timePorcion)/CLOCKS_PER_SEC << endl;

			f.close();
			cout << endl;

    }
		fin = clock();
		timetestR += fin - inicio;
		timetestC += inicio_trozos;
		//timetestR = timetestR - t_salidas;


		/*cout << "Resultado Inferencia Recursiva " << E_Par_Test.N_Examples() << endl;
		for (int t=0; t<DisparosRecursiva.size(); t++){
			cout << "(" <<t<<")\t"<<DisparosRecursiva[t].first<<"("<< DisparosRecursiva[t].second.first<<") --> " << DisparosRecursiva[t].second.second<<endl;
		}*/

		//t_salidas = 0;
    //inicio = clock();
		// Inferencia clasica    Si se descomenta (A) es suficiente con poner un /* en la siguiente linea
		// =========================================== INFERENCIA CLASICA =================================================
		//Patrones2.TestearPatronesBasicosClassic(E_Par_Test, V2, result[par]);
		//vector<pair<int,pair <string, double> > >DisparosClasica;
		//Patrones2.TestearPatronesBasicosClassicDisparos(E_Par_Test, V2, result[par], DisparosClasica);

		//fin = clock();
		//timetestC += fin - inicio;
		//timetestC = timetestC - t_salidas;

		//cout << "Resultado Inferencia Clasica " << E_Par_Test.N_Examples() << endl;
		//for (int t=0; t<DisparosClasica.size(); t++){
		//	cout << "(" <<t<<")\t"<<DisparosClasica[t].first<<"("<< DisparosClasica[t].second.first<<") --> " << DisparosClasica[t].second.second<<endl;
		//}


    /*cout << "Comparando las dos inferencias" << endl;
		double igual_clase_i=0;
		double igual_adaptacion_i=0;
		double igual_regla_i=0;
		for (int t=0; t<E_Par_Test.N_Examples(); t++){
			if (DisparosRecursiva[t].first == DisparosClasica[t].first){
				igual_clase_i++;
			}
			if (DisparosRecursiva[t].second.first == DisparosClasica[t].second.first){
				igual_regla_i++;
				//cout << "(" <<t<<")\t"<<"("<< DisparosClasica[t].second.first<<") ("<< DisparosRecursiva[t].second.first<<")"<<endl;
			}
			if (fabs(DisparosRecursiva[t].second.second - DisparosClasica[t].second.second)<0.0000001){
				//cout << "(" <<t<<")\t"<<DisparosClasica[t].first<<" | "<< DisparosRecursiva[t].first << endl;
				igual_adaptacion_i++;
			}
			//if (fabs(DisparosRecursiva[t].second.second - DisparosClasica[t].second.second)>=0.0000001){
			//	cout << "(" <<t<<")\t"<<DisparosClasica[t].second.second<<" | "<< DisparosRecursiva[t].second.second << endl;
			//}
		}
		cout << "\tMisma Clase: " << (igual_clase_i*100/E_Par_Test.N_Examples()) << "\%  Misma Regla: " << (igual_regla_i*100/E_Par_Test.N_Examples()) << "\%  Misma Adapt: " << (igual_adaptacion_i*100/E_Par_Test.N_Examples()) << "\%" <<endl;
    igual_clase += igual_clase_i*100/E_Par_Test.N_Examples();
		igual_adaptacion += igual_adaptacion_i*100/E_Par_Test.N_Examples();
		igual_regla += igual_regla_i*100/E_Par_Test.N_Examples();


    {cout << "Pulsa un tecla para continuar\n";char ch; cin >> ch;}
    */

		ProcesarResultados(resultRecur[par]);
		PintaResultadosTest(resultRecur[par], true);
		/*cerr << "% Acierto Test (Global) : " << result[par].acierto_global << endl;
		cerr << "% nuevos patrones Test  : " << result[par].porcentaje_nuevos_patrones << endl;
		cerr << "% Test sin nocubiertos  : " << result[par].acierto_sinNoCubiertos << endl;
		cerr << "%Error intrinseco       : " << result[par].error_intrinseco << endl;
		cerr << endl;*/
		//Pausa();
}

	/*cerr << endl;

	cerr << "Numero medio de Ejemplos: " << SumReSult[11] << endl;
	cerr << "Numero medio de Patrones: " << SumReSult[0] << endl;
	cerr << "         Mayoria Clase 0: " << SumReSult[6] << endl;
	cerr << "                  Empate: " << SumReSult[8] << endl;
	cerr << "         Mayoria Clase 1: " << SumReSult[7] << endl;
	cerr << "Maxima Redundancia      : " << SumReSult[9] << endl;
	cerr << "Media ejemplos/patron   : " << SumReSult[11]/SumReSult[0] << endl;
	cerr << "Mediana Redundancia     : " << SumReSult[10] << endl;
	cerr << "Patrones con 1 ejemplo  : " << SumReSult[3] << endl;
	cerr << "Patrones de 1 sola clase: " << SumReSult[1] << endl;
	cerr << "Media Train sólo 1 clase: " << SumReSult[4] << endl;
	cerr << "Media Train clase mayor : " << SumReSult[5] << endl;
	cerr << "% Acierto Test (Global) : " << SumReSult[13] << endl;
	cerr << "% nuevos patrones Test  : " << SumReSult[12] << endl;
	cerr << "% Test sin nocubiertos  : " << SumReSult[14] << endl;
	cerr << "%Error intrinseco       : " << 100 -SumReSult[14] << endl;
	cerr << endl;
*/

cerr << "\n\t------------------ Media Test -----------" << endl;
TestResult aux2;
CalcularMediaTestResult(resultRecur, aux2);
PintaResultadosTest(aux2, true);

s.clear();
s.push_back(E.N_Examples());
s.push_back(media_patrones/10.0);
s.push_back(porDefecto.acierto_global);
s.push_back(porDefecto.error_intrinseco);
s.push_back(aux2.acierto_global);
s.push_back(aux2.acierto_sinNoCubiertos);
s.push_back(aux2.error_intrinseco);
s.push_back(aux2.porcentaje_nuevos_patrones);
s.push_back( (1.0*timetrain/num_par) / CLOCKS_PER_SEC );
s.push_back( (1.0*timetestR/num_par) / CLOCKS_PER_SEC ); // Inferencia no clásica
s.push_back( (1.0*timetestC/num_par) / CLOCKS_PER_SEC ); // Inferencia Clásica
s.push_back(igual_clase/num_par);  // Porcentaje en el que se disparan reglas de devuelven la misma clase en las dos inferencias
s.push_back(igual_regla/num_par);  // Porcentaje en el que se dispara la misma regla en las dos inferencias
s.push_back(igual_adaptacion/num_par);  // Porcentaje en el que coincide la adaptacion de las dos regla en las dos inferencias

salida[0] = s;

for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
	salida[i+1].push_back(100.0 * aux2.acc[i]/(aux2.cubiertos[i]+aux2.no_cubiertos[i]));
	if (aux2.cubiertos[i]>0){
		salida[i+1].push_back(100.0 * aux2.acc[i]/(aux2.cubiertos[i]));
		salida[i+1].push_back(100.0 - 100.0 * aux2.acc[i]/(aux2.cubiertos[i]));
	}
	else {
		salida[i+1].push_back(0);
		salida[i+1].push_back(100);
	}
	salida[i+1].push_back(100.0 * aux2.no_cubiertos[i]/(aux2.cubiertos[i]+aux2.no_cubiertos[i]));
}


//Pausa();

//exit(0);

return salida;

}




//-----------------------------------------------------------------------------------------------------

vector<vector<double> > CalculoEstimadoSobrePatrones3(example_set & E, const VectorVar &V, int num_par){
  vector<TestResult> result;
	TestResult porDefecto, aux, porDefecto2;
	VectorVar V2, V3;
	V2 = V;
	V3 = V;

	// Creo el dominio con 2 etiquetas en las variables continuas
	for (int i = 0; i<V.N_Antecedente(); i++){
		if (!V.Variable(i).IsDiscrete()){
			variable_t aux;
			aux.Asigna(2, V.Variable(i).Inf_Range(), V.Variable(i).Sup_Range(), true, true, V.Variable(i).Name());
			V2.Asigna(i, aux);
		}
		else {
			V2.Asigna(i, V.Variable(i));
		}
	}

	// Creo el dominio con 3 etiquetas en las variables continuas
	for (int i = 0; i<V.N_Antecedente(); i++){
		if (!V.Variable(i).IsDiscrete()){
			variable_t aux;
			aux.Asigna(3, V.Variable(i).Inf_Range(), V.Variable(i).Sup_Range(), true, true, V.Variable(i).Name());
			V3.Asigna(i, aux);
		}
		else {
			V3.Asigna(i, V.Variable(i));
		}
	}



	Pattern P2(V2), P3(V3);

	vector <vector <double> > salida;
	vector<double> s;

  s.clear();
	for (int i=0; i<V.SizeDomain(V.Consecuente())+1; i++){
		salida.push_back(s);
	}

	P2.ExtraerPatronesBasicos(E,V2, porDefecto);
	P3.ExtraerPatronesBasicos(E,V3, porDefecto2);
	cerr << "Numero de Ejemplos: " << E.N_Examples() << endl;
	for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
		cerr << "           clase " << i << ": " << P2.N_Ejemplos_paraClase(i) << " (" <<100.0*P2.N_Ejemplos_paraClase(i) / E.N_Examples() << ")" <<endl;
		cerr << "           clase " << i << ": " << P3.N_Ejemplos_paraClase(i) << " (" <<100.0*P3.N_Ejemplos_paraClase(i) / E.N_Examples() << ")" <<endl;
		salida[i+1].push_back(P2.N_Ejemplos_paraClase(i)+P3.N_Ejemplos_paraClase(i));
	}
	cerr << "Numero de Patrones: " << P2.N_Pattern()+P3.N_Pattern() << endl;
	for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
		cerr << "           clase " << i << ": " << P2.N_Patrones_paraClase(i) << " (" <<100.0*P2.N_Patrones_paraClase(i) / P2.N_Pattern() << ")" <<endl;
		cerr << "           clase " << i << ": " << P3.N_Patrones_paraClase(i) << " (" <<100.0*P3.N_Patrones_paraClase(i) / P3.N_Pattern() << ")" <<endl;
		salida[i+1].push_back(P2.N_Patrones_paraClase(i)+P2.N_Patrones_paraClase(i));
	}
	cerr << "Media ejemplos/patron(2)   : " << 1.0*E.N_Examples()/P2.N_Pattern() << endl;
	cerr << "Media ejemplos/patron(3)   : " << 1.0*E.N_Examples()/P3.N_Pattern() << endl;
	//P.TestearPatronesBasicos(E,V,porDefecto);
	PintaResultadosTraining(porDefecto);

	for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
		salida[i+1].push_back(100.0 * porDefecto.acc[i]/(porDefecto.cubiertos[i]));
		salida[i+1].push_back(100.0 - 100.0 * porDefecto.acc[i]/(porDefecto.cubiertos[i]));
	}

	//Pausa();


	for (int par=0; par<num_par; par++){
		Pattern Patrones2(V2);
		Pattern Patrones3(V3);
		cout <<"\n===============\n Ejecucion " << par << "\n===============\n";

		example_set E_Par_Completo=E.Extract_Training_Set(par, V.Consecuente(), 0);
		example_set E_Par_Test = E.Extract_Test_Set(par, V.Consecuente());

		result.push_back(porDefecto);
		Patrones2.ExtraerPatronesBasicos(E_Par_Completo,V2, result[par]);
		//Patrones3.ExtraerPatronesBasicos(E_Par_Completo,V3, result[par]);
		//Patrones.TestearPatronesBasicos(E_Par_Completo,V,result[par]);

		cerr << "Numero de Ejemplos: " << E_Par_Completo.N_Examples() << endl;
		for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
			cerr << "         clase " << i << ": " << Patrones2.N_Ejemplos_paraClase(i) << " (" <<100.0*Patrones2.N_Ejemplos_paraClase(i) / E_Par_Completo.N_Examples() << ")" <<endl;
			//cerr << "         clase " << i << ": " << Patrones3.N_Ejemplos_paraClase(i) << " (" <<100.0*Patrones3.N_Ejemplos_paraClase(i) / E_Par_Completo.N_Examples() << ")" <<endl;
		}
		cerr << "Numero de Patrones: " << Patrones2.N_Pattern()+Patrones3.N_Pattern() << endl;
		for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
			cerr << "         clase " << i << ": " << Patrones2.N_Patrones_paraClase(i) << " (" <<100.0*Patrones2.N_Patrones_paraClase(i) / Patrones2.N_Pattern() << ")" <<endl;
			//cerr << "         clase " << i << ": " << Patrones3.N_Patrones_paraClase(i) << " (" <<100.0*Patrones3.N_Patrones_paraClase(i) / Patrones3.N_Pattern() << ")" <<endl;
		}
		cerr << "Media ejemplos/patron   : " << 1.0*E_Par_Completo.N_Examples()/(Patrones2.N_Pattern()+Patrones3.N_Pattern()) << endl;
		for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
			cerr << "         clase " << i << ": " << 1.0 * Patrones2.N_Ejemplos_paraClase(i) / Patrones2.N_Patrones_paraClase(i) <<endl;
			//cerr << "         clase " << i << ": " << 1.0 * Patrones3.N_Ejemplos_paraClase(i) / Patrones3.N_Patrones_paraClase(i) <<endl;
		}
		cerr << "Porcentaje sobre el global de patrones: " << 100.0 * (Patrones2.N_Pattern()+Patrones3.N_Pattern()) / (P2.N_Pattern()+P3.N_Pattern()) << endl;
		for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
			cerr << "         clase " << i << ": " << Patrones2.N_Patrones_paraClase(i) << " (" <<100.0*Patrones2.N_Patrones_paraClase(i) / (P2.N_Pattern()+P3.N_Pattern()) << ")" <<endl;
			//cerr << "         clase " << i << ": " << Patrones3.N_Patrones_paraClase(i) << " (" <<100.0*Patrones3.N_Patrones_paraClase(i) / (P2.N_Pattern()+P3.N_Pattern()) << ")" <<endl;
		}
		//cerr << "Media Train clase mayor : " << result[par].acierto_global << endl;
		cerr << "\t-------------------- Training -----------" << endl;
		PintaResultadosTraining(result[par]);

		cerr << "\n\t--------------------- Test  -----------" << endl;
		double peso1, peso2;
		result[par].cubiertos.clear();
	  result[par].no_cubiertos.clear();
	  result[par].acc.clear();
	  result[par].error_intrinseco_porClase.clear();

	  for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
	    result[par].cubiertos.push_back(0);
	    result[par].no_cubiertos.push_back(0);
	    result[par].acc.push_back(0);
	    result[par].error_intrinseco_porClase.push_back(0);
	  }

		for (int i=0; i<E_Par_Test.N_Examples(); i++){
			int cl1 = Patrones2.TestearUnEjemploPatronesBasicos(E_Par_Test, i, V2, peso1);
			//int cl2 = Patrones3.TestearUnEjemploPatronesBasicos(E_Par_Test, i, V3, peso2);
			int cl2 = cl1;
			int cl;

			// Toma de decisiones
			if (cl1 == cl2){
				cl = cl1;
			}
			else if (cl1 == -1){
				cl = cl2;
			}
			else if (cl2 == -1){
				cl = cl1;
			}
			else{
				if (peso2 >= peso1){
					cl = cl2;
				}
				else {
					cl = cl1;
				}
			}

			int clase = E_Par_Test.Data(i,V.Consecuente());
			if (cl == -1)	{
				result[par].no_cubiertos[clase]++;
			}
			else {
				result[par].cubiertos[clase]++;
				if (cl == clase){
					result[par].acc[clase]++;
				}
			}

		}
		PintaResultadosTest(result[par], true);
		/*cerr << "% Acierto Test (Global) : " << result[par].acierto_global << endl;
		cerr << "% nuevos patrones Test  : " << result[par].porcentaje_nuevos_patrones << endl;
		cerr << "% Test sin nocubiertos  : " << result[par].acierto_sinNoCubiertos << endl;
		cerr << "%Error intrinseco       : " << result[par].error_intrinseco << endl;*/
		cerr << endl;
		//Pausa();
}

	/*cerr << endl;

	cerr << "Numero medio de Ejemplos: " << SumReSult[11] << endl;
	cerr << "Numero medio de Patrones: " << SumReSult[0] << endl;
	cerr << "         Mayoria Clase 0: " << SumReSult[6] << endl;
	cerr << "                  Empate: " << SumReSult[8] << endl;
	cerr << "         Mayoria Clase 1: " << SumReSult[7] << endl;
	cerr << "Maxima Redundancia      : " << SumReSult[9] << endl;
	cerr << "Media ejemplos/patron   : " << SumReSult[11]/SumReSult[0] << endl;
	cerr << "Mediana Redundancia     : " << SumReSult[10] << endl;
	cerr << "Patrones con 1 ejemplo  : " << SumReSult[3] << endl;
	cerr << "Patrones de 1 sola clase: " << SumReSult[1] << endl;
	cerr << "Media Train sólo 1 clase: " << SumReSult[4] << endl;
	cerr << "Media Train clase mayor : " << SumReSult[5] << endl;
	cerr << "% Acierto Test (Global) : " << SumReSult[13] << endl;
	cerr << "% nuevos patrones Test  : " << SumReSult[12] << endl;
	cerr << "% Test sin nocubiertos  : " << SumReSult[14] << endl;
	cerr << "%Error intrinseco       : " << 100 -SumReSult[14] << endl;
	cerr << endl;
*/

cerr << "\n\t------------------ Media Test -----------" << endl;
TestResult aux2;
CalcularMediaTestResult(result, aux2);
PintaResultadosTest(aux2, true);

s.clear();
s.push_back(E.N_Examples());
s.push_back(P2.N_Pattern()+P2.N_Pattern());
s.push_back(porDefecto.acierto_global);
s.push_back(porDefecto.error_intrinseco);
s.push_back(aux2.acierto_global);
s.push_back(aux2.acierto_sinNoCubiertos);
s.push_back(aux2.error_intrinseco);
s.push_back(aux2.porcentaje_nuevos_patrones);
salida[0] = s;

for (int i=0; i<V.SizeDomain(V.Consecuente()); i++){
	salida[i+1].push_back(100.0 * aux2.acc[i]/(aux2.cubiertos[i]+aux2.no_cubiertos[i]));
	if (aux2.cubiertos[i]>0){
		salida[i+1].push_back(100.0 * aux2.acc[i]/(aux2.cubiertos[i]));
		salida[i+1].push_back(100.0 - 100.0 * aux2.acc[i]/(aux2.cubiertos[i]));
	}
	else {
		salida[i+1].push_back(0);
		salida[i+1].push_back(100);
	}
	salida[i+1].push_back(100.0 * aux2.no_cubiertos[i]/(aux2.cubiertos[i]+aux2.no_cubiertos[i]));
}


//Pausa();

//exit(0);

return salida;

}







int main(int argc, char *argv[])
{

	for (char i=0; i<10; i++){
		simbolo.push_back('0'+i);
	}
	for (char i=0; i<25; i++){
		simbolo.push_back('A'+i);
	}
	for (char i=0; i<25; i++){
		simbolo.push_back('a'+i);
	}
	for (char i=0; i<40; i++){
		simbolo.push_back(200+i);
	}


	if (argc < 2)
	{
		cout << "Number of parameters incorrect....\n";
		cout << "Sintaxis:\n";
		cout << "aprendizajegenetico <semilla_problema>\n";
		cout << "       or\n";
		cout << "aprendizajegenetico <semilla_problema> implicity\n\n";
		exit(1);
	}

	// Numero de particiones
	int num_par=10;

	string fichname, fich_extension, comando;

	fichname = argv[1];
	fich_extension = fichname;
	fich_extension = fich_extension + ".dom";
	VectorVar V(fich_extension.c_str());
	example_set E;
	int seed = 0;
	int nLab = 2;

	if (argc==2)
	{
		fich_extension = fichname;
		fich_extension = fich_extension + ".datos";
		E.AddExampleFich(fich_extension.c_str(),1,1);
	}
	else if (argc>=3)
	{
		seed = atoi(argv[2]);
		int item_par=0;
		fich_extension = fichname;
		comando="0";
		comando[0]=comando[0]+item_par;
		fich_extension = fich_extension+comando;
		fich_extension = fich_extension+".datos";
		int x = 0;
		if (argc > 3) {
			nLab = atoi(argv[3]);
		}
		cout << "Cargando el fichero de datos ..." << endl;
		while (Existe(fich_extension))
		{
			E.AddExampleFich(fich_extension.c_str(),item_par,x); // Ultimo componente false=orden original, true=orden aleatorio
			cout << "...... cargado " << fich_extension << endl;
			item_par++;
			fich_extension = fichname;
			comando="0";
			comando[0]=comando[0]+item_par;
			fich_extension = fich_extension+comando;
			fich_extension = fich_extension+".datos";
		}
		num_par=item_par;
		//E.OrderByClass(1, V.Consecuente());
	}

	cout << "Fin del proceso de carga ..." << endl;
	fichname = extract_noun(fichname);

	//int nn = rand() % 1000;
	fich_extension = fichname + itos(seed) + ".res";
	ofstream out(fich_extension.c_str());
	out << "Experiment: " << fichname << endl;
	out.close();

	comando = "date >> ";
	comando = comando+fich_extension;
	system (comando.c_str());

	out.open(fich_extension.c_str(),ios::app);
	out << "===================================================\n\n";
	out.close();

	srand(seed);
	cout << "SEED " << seed << endl;



	vectordouble a1, a2(5), a3;


	//R.Load(nomfich);
	//numero de poblaciones binarias
	int NPobBin=2;
	//numero de poblaciones enteras
	int NPobEnt=2;
	//numero de poblaciones reales
	int NPobRea=2;               //5
	//Elite de la poblacion genetica
	int Elite=10;
	//numero de componentes de la funcion fitness
	int NFitness=7;
	//numero de individuos de la poblacion
	int Indiv=100;
	//numero total de poblaciones del genetico
	int NPobTotal=NPobBin+NPobEnt+NPobRea;
	//Establece el rango de variacion de cada poblacion entera [0,rango[i]]
	int *rango=new int[NPobEnt];
	//Establece el rango de variacion de cada poblacion real [rangoi[i],rangos[i]]
	double *rangoi=new double[NPobRea];
	double *rangos=new double[NPobRea];
	//Establece los parametros del genetico:
	//  Tamano de cada poblacion, prob. de mut y prob. de cruce.
	int *tama=new int[NPobTotal];
	double *mut=new double[NPobTotal];
	double *cruce =new double[NPobTotal];

	//Configuracion del Genetico
	for (int i=0; i<NPobBin; i++)
	{
		mut[i]=0.01;   // 0.05
		cruce[i]=0.3; // 0.6
	}

	for (int i=NPobBin; i<NPobBin+NPobEnt; i++)
	{
		mut[i]=0.01;
		cruce[i]=0.0;
	}

	mut[NPobBin]=0.01;  // 0.01
	cruce[NPobBin]=0.2;


	for (int i=NPobBin+NPobEnt; i<NPobBin+NPobEnt+NPobRea; i++)
	{
		mut[i]=0.05;     // 0.01
		cruce[i]=0;
	}


	// Valores para el nivel de variable
	mut[NPobBin+NPobEnt]=0.5/V.N_Antecedente();     // 0.01
	cruce[NPobBin+NPobEnt]=0.2;
	mut[1]=0.1;


	int antecedente;

	string seedfich=fichname;
	string ext=".rl0";
	string ext3=".rl_descrip0";
	string ext2=".rules";
	string nomfich, nomfich2, nomfich3;

	if (argc==2)
	E.Generate_Partitions(num_par);


	example_set E_Par, E_Par_Test, E_Par_Completo;

	//==================================================
  // Calculo de Patrones
	//==================================================
  for (int numLabel = 3; numLabel<4; numLabel++){
	  int num_variables_continuas = 0;
		for (int i = 0; i<V.N_Antecedente(); i++){
			if (!V.Variable(i).IsDiscrete()){
				num_variables_continuas++;
			}
		}
		if (false /*num_variables_continuas == 0*/) {
			vector<vector<double> > salida;
			salida = CalculoEstimadoSobrePatrones2_Para5subconjuntos(E, V, V, V, num_par);

			fstream f;
			f.open("patrBasicos.csv", std::fstream::out | std::fstream::app);

			for (int j=0; j<salida.size(); j++){
				f << seedfich << " 0 0";
				if (j==0){
					f << " Global";
				}
				else {
					f << " Class" << j;
				}
				for (int i=0; i<salida[j].size(); i++){
					f << " " << salida[j][i];
				}
				f << endl;
			}
			f.close();
		}
		else {
			// Formato 2-3
			fstream f;
			vector<vector<double> > salida;
			//salida = CalculoEstimadoSobrePatrones3(E, V,  num_par);

	    int numLabel = nLab;
			salida = CalculoEstimadoSobrePatronesConNEtiquetas(numLabel, E, V, num_par);


			for (int j=0; j<salida.size(); j++){
				cout << seedfich << " " << numLabel << " " << num_variables_continuas;
				if (j==0){
					cout << " Global";
				}
				else {
					cout << " Class" << j;
				}
				for (int i=0; i<salida[j].size(); i++){
					cout << " " << salida[j][i];
				}
				cout << endl;
			}


			f.open("patrBasicos.csv", std::fstream::out | std::fstream::app);

			for (int j=0; j<salida.size(); j++){
				f << seedfich << " " << numLabel << " " << num_variables_continuas;
				if (j==0){
					f << " Global";
				}
				else {
					f << " Class" << j;
				}
				for (int i=0; i<salida[j].size(); i++){
					f << " " << salida[j][i];
				}
				f << endl;
			}
			f.close();
		}

		/*for (int k=5; k<6; k++){
			VectorVar V2 = V;
			VectorVar V3 = V;
			VectorVar V5 = V;


			V5.Asigna(V.Consecuente(), V.Variable(V.Consecuente()));

			for (int i = 0; i<V.N_Antecedente(); i++){
				if (!V.Variable(i).IsDiscrete()){
					variable_t aux;
					aux.Asigna(5, V.Variable(i).Inf_Range(), V.Variable(i).Sup_Range(), true, true, V.Variable(i).Name());
					V5.Asigna(i, aux);
				}
				else {
					V5.Asigna(i, V.Variable(i));
				}
			}


			for (int i = 0; i<V.N_Antecedente(); i++){
				if (!V.Variable(i).IsDiscrete()){
					variable_t aux;
					aux.Asigna(3, V.Variable(i).Inf_Range(), V.Variable(i).Sup_Range(), true, true, V.Variable(i).Name());
					V3.Asigna(i, aux);
				}
				else {
					V3.Asigna(i, V.Variable(i));
				}
			}


			for (int i = 0; i<V.N_Antecedente(); i++){
				if (!V.Variable(i).IsDiscrete()){
					variable_t aux;
					aux.Asigna(2, V.Variable(i).Inf_Range(), V.Variable(i).Sup_Range(), true, true, V.Variable(i).Name());
					V2.Asigna(i, aux);
				}
				else {
					V2.Asigna(i, V.Variable(i));
				}
			}

			//V2.Pinta();
			//Pausa();

			//salida = CalculoEstimadoSobrePatrones2_Para5subconjuntos(E, V5, V3, V2, num_par);

			salida = CalculoEstimadoSobrePatrones3(E, V,  num_par);


			f.open("patrBasicos.csv", std::fstream::out | std::fstream::app);

			for (int j=0; j<salida.size(); j++){
				f << seedfich << " " << k << " " << num_variables_continuas;
				if (j==0){
					f << " Global";
				}
				else {
					f << " Class" << j;
				}
				for (int i=0; i<salida[j].size(); i++){
					f << " " << salida[j][i];
				}
				f << endl;
			}
			f.close();
		}*/
	}

  delete [] rango;
	delete [] rangoi;
	delete [] rangos;
	delete [] tama;
	delete [] mut;
	delete [] cruce;
}



bool vectoresIguales_E(double *v1, double *v2, int tama){
  int i=0;
/*  for (int t=0; t<tama; t++){
    cout << v1[t] << " | " << v2[t] << endl;
  }*/

  while (i<tama and (v1[i] == v2[i] /*or (i==1 and abs(v1[i]-v2[i])<=1)*/ ) ){
    i++;
  }

  return (i==tama);
}

bool vectoresIguales(double *v1, double *v2, int tama){
  int i=0;
  /*for (int t=0; t<tama; t++){
    cout << v1[t] << " | " << v2[t] << endl;
  }*/

  while (i<tama and v1[i]==v2[i] ){
    i++;
  }

  return (i==tama);
}


int DistanciaHamming(const string &s1, const string &s2){
  int cont=0;
  for (int i=0; i<s1.length(); i++)
    if (s1[i]!=s2[i] and s1[i]!='?' and s2[i]!='?')
      cont++;

  return cont;
}


// Devuelve la primera diferencia entre 2 patrones a partir de pos_inicial
int DiferenciaDistanciaHamming(const string &s1, const string &s2, int pos_inicial){
  int i=0;
  for (int i=0; i<s1.length(); i++)
  while (i <s1.length() and  s1[i]!=s2[i] and s1[i]!='?' and s2[i]!='?')
      i++;

  if (i==s1.length())
    return -1;
  else
    return i;
}


double Evento1_Especificar(multipoblacion &G, const VectorVar &V, example_set &E_Par, ruleset &R, ruleset &R_dis, ruleset &R_mix, int *ponderacion_por_clase, VectorVar &F,
       double ***Tabla, double ***Tabla_Nuevas, int i, double *multi_valor){
  double valor1[10];
  int k=G.N_individuos()-1;
  genetcode *aux, *code_dis, *code_des;
  int NFitness=G.N_Valoracion();
  double peso_anterior1, peso_anterior2;
  double nmas1, nmas2, nmenos1, nmenos2, nmas, nmenos;

  bool discriminante=false, reducir_reglas=true;

  //for (int i=0; i<R.N_rule(); i++){

    //cout << "Patron antes del proceso:\n";
    //string patr1 = Patron2Reglas(R,R_dis,V,i);

    aux = R_mix.Get_Code_Rule(i);
    int tam = aux->TamaBinary(0);
    char *cad= new char[tam+1];
    for (int ss=0; ss<tam; ss++)
      cad[ss]= aux->GetValueBinary(0,ss);
    cad[tam]='\0';
    string patr1 = cad;
    delete [] cad;
    //cout << patr1 << endl;

    // saco la regla descriptiva
    code_des=R.Extract(i, peso_anterior2, nmas2, nmenos2);

    // copio de regla discriminante
    nmas = R_dis.Get_N_Mas(i);
    nmenos= R_dis.Get_N_Menos(i);
    aux = R_dis.Get_Code_Rule(i);
    *code_des = *aux;

    // Aniado la regla discriminante a la poblacion en la posicion k
    G.PutCodigo(k,*code_des);

    //R.Insert(i,code_des, nmas2, nmenos2, false);

    genetcode* aux_R =R_dis.Get_Code_Rule(i);
    int len = aux_R->TamaBinary(0);
    cad = new char[len+1];
    for (int j=0; j<len; j++)
      cad[j]=aux_R->GetValueBinary(0,j);
    cad[len]='\0';
    string regla = cad;
    delete [] cad;

    //cout << "Regla inicial: \n" << regla << endl;


    //cout << "\tAcierto antes de Especificar: " <<Aciertos(R,V,E_Par, F) << endl;

    int *basica_para_regla = new int[R.N_rule()+1];
    int *errores_para_regla = new int[R.N_rule()+1];
    Obtener_Aciertos_y_Errores_Reglas_Aprendidas(E_Par, basica_para_regla, errores_para_regla, R.N_rule()+1);
    int *reglas_a_eliminar = new int[R.N_rule()+1];
    int n_reglas_a_eliminar;
    double peso1 = peso_anterior1;

    Bondad_E(G,V,E_Par,R,valor1, NFitness, k, discriminante, reducir_reglas, ponderacion_por_clase,
		      R.N_rule(), Tabla, Tabla_Nuevas, F, basica_para_regla, errores_para_regla, reglas_a_eliminar,
		      n_reglas_a_eliminar, true, regla, nmas, nmenos);

    Weight_Individual(G,V,E_Par,F,Tabla, Tabla_Nuevas, k, nmas1, nmenos1);

    //cout << "Especificar\n";
    //cout << patr1 << endl;
    //cout << "Antes:   " << valor1[0] << ", " << valor1[1] << ", "<< valor1[2] << ", " << valor1[3] << endl;


    Obtener_Regla_Descriptiva_E(G,V,E_Par,R,k,discriminante, reducir_reglas,ponderacion_por_clase,Tabla,
				Tabla_Nuevas, F, basica_para_regla, errores_para_regla,
				reglas_a_eliminar, n_reglas_a_eliminar, false, patr1, false, regla, nmas, nmenos);


    peso1 =Bondad_E(G,V,E_Par,R,valor1, NFitness, k, discriminante, reducir_reglas, ponderacion_por_clase,
		      R.N_rule(), Tabla, Tabla_Nuevas, F, basica_para_regla, errores_para_regla, reglas_a_eliminar,
		      n_reglas_a_eliminar,true, regla, nmas, nmenos);

    peso1 = Weight_Individual(G,V,E_Par,F,Tabla, Tabla_Nuevas, k, nmas2, nmenos2);

    //cout << "Despues: " << valor1[0] << ", " << valor1[1] << ", " << valor1[2] << ", " << valor1[3] << " Peso= " << peso1 << endl << endl;


    //code_des=R.Extract(i, peso_anterior1, nmas1, nmenos1);
    G.Codigo(k,*code_des);
    R.Insert(i,code_des,nmas2,nmenos2, false);
    //R.Recalculate_Weight(i);

    //cout << "Peso de la regla: " << R.Get_Weight(i);

      delete [] basica_para_regla;
      delete [] errores_para_regla;
      delete [] reglas_a_eliminar;

    //cout << "Patron despues del proceso:\n";
    //string patr2 = Patron2Reglas(R,R_dis,V,i);
    //cout << patr2 << endl;
    //Comparar2patrones(patr1,patr2);
    //cout << endl;

    aux_R =R.Get_Code_Rule(i);
    len = aux_R->TamaBinary(0);
    cad = new char[len+1];
    for (int j=0; j<len; j++)
      cad[j]=aux_R->GetValueBinary(0,j);
    cad[len]='\0';
    regla = cad;
    delete [] cad;

    //cout << "Regla resultante: \n" << regla << endl;
    //string patr2 = Patron2Reglas(R,R_dis,V,i);
    /*aux = R_mix.Get_Code_Rule(i);
    for (int j=0; j<aux->TamaBinary(0); j++)
      aux->PutValueBinary(0,j,patr2[j]);*/

  //}
  //char ch;
  //cin >> ch;

    multi_valor[0]=valor1[0];
    multi_valor[1]=valor1[1];
    multi_valor[2]=valor1[2];
    multi_valor[3]=0;
    //multi_valor[3]=valor1[3];

    //cout << "\tAcierto despues de Especificar: " << Aciertos(R,V,E_Par, F) << endl;

  return valor1[0];
}



double Evento1_Generalizar(multipoblacion &G, const VectorVar &V, example_set &E_Par, ruleset &R, ruleset &R_dis, ruleset &R_mix, int *ponderacion_por_clase, VectorVar &F,
       double ***Tabla, double ***Tabla_Nuevas, int i, double *multi_valor){
  double peso_anterior1, peso_anterior2;
  double valor1[10];
  int k=G.N_individuos()-1;
  genetcode *aux, *code_dis, code_des;
  int NFitness=G.N_Valoracion();
  double nmas1, nmas2, nmenos1, nmenos2, nmas, nmenos;

  bool discriminante=false, reducir_reglas=true;

  //for (int i=0; i<R.N_rule(); i++){

    //cout << "Patron antes del proceso:\n";
    //string patr1 = Patron2Reglas(R,R_dis,V,i);
    //cout << patr1 << endl;

    genetcode* aux_mix =R_mix.Get_Code_Rule(i);
    int len = aux_mix->TamaBinary(0);
    char *cad = new char[len+1];
    for (int j=0; j<len; j++)
      cad[j]=aux_mix->GetValueBinary(0,j);
    cad[len]='\0';
    string patr1 = cad;
    delete [] cad;

    // Saco la regla discriminante de la base de reglas descriptivas y la elimino
    code_dis=R_dis.Extract(i, peso_anterior2);

    // saco la regla descriptiva
    //code_des=R.Extract(i, peso_anterior2, nmas2, nmenos2);

    peso_anterior1 = R.Get_Weight(i);
    nmas = R.Get_N_Mas(i);
    nmenos = R.Get_N_Menos(i);
    aux = R.Get_Code_Rule(i);
    code_des = *aux;

    // Aniado la regla descriptiva a la poblacion en la posicion k
    G.PutCodigo(k,code_des);


    genetcode* aux_R =R.Get_Code_Rule(i);
    len = aux_R->TamaBinary(0);
    cad = new char[len+1];
    for (int j=0; j<len; j++)
      cad[j]=aux_R->GetValueBinary(0,j);
    cad[len]='\0';
    string regla = cad;
    delete [] cad;

    genetcode *code_R = R.Extract(i,peso_anterior1);

    //cout << "Regla inicial: \n" << regla << endl;

    //cout << "\tAcierto antes de Generalizar: " << Aciertos(R,V,E_Par, F) << endl;

    int *basica_para_regla = new int[R.N_rule()+1];
    int *errores_para_regla = new int[R.N_rule()+1];
    Obtener_Aciertos_y_Errores_Reglas_Aprendidas(E_Par, basica_para_regla, errores_para_regla, R.N_rule()+1);
    int *reglas_a_eliminar = new int[R.N_rule()+1];
    int n_reglas_a_eliminar;

    Bondad_E(G,V,E_Par,R,valor1, NFitness, k, discriminante, reducir_reglas, ponderacion_por_clase,
		      R.N_rule(), Tabla, Tabla_Nuevas, F, basica_para_regla, errores_para_regla, reglas_a_eliminar,
		      n_reglas_a_eliminar, true, regla, nmas, nmenos);


    Weight_Individual(G,V,E_Par,F,Tabla, Tabla_Nuevas, k, nmas1, nmenos1);

    //cout << "Generalizar\n";
    //cout << patr1 << endl;
    //cout << "Antes:   " << valor1[0] << ", " << valor1[1] << ", " << valor1[2] << ", " << valor1[3] << endl;


    Obtener_Regla_Discriminante_E(G,V,E_Par,R,k,discriminante, reducir_reglas,ponderacion_por_clase,Tabla,
				Tabla_Nuevas, F, basica_para_regla, errores_para_regla,
				reglas_a_eliminar, n_reglas_a_eliminar, false, patr1, false, regla, nmas, nmenos);


    double peso2 =Bondad_E(G,V,E_Par,R,valor1, NFitness, k, discriminante, reducir_reglas, ponderacion_por_clase,
		      R.N_rule(), Tabla, Tabla_Nuevas, F, basica_para_regla, errores_para_regla, reglas_a_eliminar,
		      n_reglas_a_eliminar, true, regla, nmas, nmenos);


    peso2 = Weight_Individual(G,V,E_Par,F,Tabla, Tabla_Nuevas, k, nmas2, nmenos2);

    //cout << "Despues: " << valor1[0] << ", " << valor1[1] << ", " << valor1[2] << ", " << valor1[3] << endl << endl;

      G.Codigo(k,*code_dis);
      R_dis.Insert(i,code_dis, nmas2, nmenos2, false);
      R_dis.Recalculate_Weight(i);

    aux_R =R_dis.Get_Code_Rule(i);
    len = aux_R->TamaBinary(0);
    cad = new char[len+1];
    for (int j=0; j<len; j++)
      cad[j]=aux_R->GetValueBinary(0,j);
    cad[len]='\0';
    regla = cad;
    delete [] cad;

    //cout << "Regla resultante: \n" << regla << endl;


      //aux=R.Extract(i, peso_anterior2);
      //*aux=code_des;
      R.Insert(i,code_R, nmas, nmenos, false);
      /*R.Anular_N_Mas_y_N_Menos(i);
      R.Add_N_Mas(i,nmas1);
      R.Add_N_Menos(i,nmenos2);
      R.Recalculate_Weight(i);      */

      delete [] basica_para_regla;
      delete [] errores_para_regla;
      delete [] reglas_a_eliminar;

    //cout << "Patron despues del proceso:\n";
    string patr2 = Patron2Reglas(R,R_dis,V,i);
    //cout << patr2 << endl;
    //Comparar2patrones(patr1,patr2);
    //cout << endl;



    //Pausa();

    /*aux = R_mix.Get_Code_Rule(i);
    for (int j=0; j<aux->TamaBinary(0); j++)
      aux->PutValueBinary(0,j,patr2[j]);*/

  //}

  //char ch;
  //cin >> ch;
    multi_valor[0]=valor1[0];
    multi_valor[1]=valor1[1];
    multi_valor[2]=valor1[2];
    multi_valor[3]=0;
    //multi_valor[3]=valor1[3];

    //cout << "\tAcierto despues de generalizar: " << Aciertos(R,V,E_Par, F) << endl;

  return valor1[0];

}



void PintaDiferenciaPatrones(const string &cad1, const string &cad2){
  int n_cambios=0;
  bool compatibles=true;
  for (int s=0; s<cad1.length(); s++)
    if (cad1[s]==cad2[s]){
      cout << "_";
    }
    else{
      if ( (cad1[s]=='x' and cad2[s]=='0') or (cad1[s]=='0' and cad2[s]=='x') ){
        cout << "X";
	compatibles=false;
      }
      else if (cad1[s]=='0' and cad2[s]=='1'){
	cout << ">";
	compatibles=false;
      }
      else if (cad1[s]=='1' and cad2[s]=='0'){
	cout << "<";
	compatibles=false;
      }
      else{
	cout << "$";
      }

      n_cambios++;
    }

  cout << endl;
}


int Comparar2patrones(const string &cad1, const string &cad2){
  int n_cambios=0;
  bool compatibles=true;
  for (int s=0; s<cad1.length(); s++)
    if (cad1[s]==cad2[s]){
      //cout << "_";
    }
    else{
      if (cad1[s]=='0' and cad2[s]=='1'){
	//cout << ">";
	compatibles=false;
      }
      else if (cad1[s]=='1' and cad2[s]=='0'){
	//cout << "<";
	compatibles=false;
      }
      else{
	//cout << "$";
      }

      n_cambios++;
    }

  //cout << endl;

  if (n_cambios==0)
    return 0; // Son iguales
  else if (compatibles)
    return 1;
  else
    return 2;
}



int Comparar2patronesPinta(const string &cad1, const string &cad2){
  int n_cambios=0;
  bool compatibles=true;
  for (int s=0; s<cad1.length(); s++)
    if (cad1[s]==cad2[s]){
      cout << "_";
    }
    else{
      if ( (cad1[s]=='x' and cad2[s]=='0') or (cad1[s]=='0' and cad2[s]=='x') ){
        cout << "X";
	compatibles=false;
      }
      else if (cad1[s]=='0' and cad2[s]=='1'){
	cout << ">";
	compatibles=false;
      }
      else if (cad1[s]=='1' and cad2[s]=='0'){
	cout << "<";
	compatibles=false;
      }
      else{
	cout << "$";
      }

      n_cambios++;
    }

  cout << endl;

  if (n_cambios==0)
    return 0; // Son iguales
  else if (compatibles)
    return 1;
  else
    return 2;
}




// Hace la union de las dos reglas y la operacion la almacena en la regla r1
void UnionReglasDescriptivas(const VectorVar &V, ruleset &R1, int r1, ruleset &R2, int r2){
  genetcode r_1, r_2, *aux;
  aux=R1.Get_Code_Rule(r1);
  r_1=*aux;
  aux=R2.Get_Code_Rule(r2);
  r_2=*aux;

 int b,*nb,n,*nn,**nnn,r,*nr;
 char **mb;
 double **mr;
 string regla1, regla2, regla3;

 r_1.GetBinary(b,nb,mb);
 r_1.GetReal(r,nr,mr);

 char *s= new char[nb[0]+1];
 for (int j=0; j<nb[0]; j++)
   s[j]=mb[0][j];
 s[nb[0]]='\0';
 regla1=s;
 delete []s;

 int longitud, tama=0;

 for (int j=0;j<V.N_Antecedente(); j++){
     longitud = V.SizeDomain(j);
     if (mr[0][j]<mr[0][nr[0]-1]){
       for (int s=0; s<longitud; s++)
	 regla1[tama+s]='1';
     }
   tama+=longitud;
 }


 r_2.GetBinary(b,nb,mb);
 r_2.GetReal(r,nr,mr);

 s= new char[nb[0]+1];
 for (int j=0; j<nb[0]; j++)
   s[j]=mb[0][j];
 s[nb[0]]='\0';
 regla2=s;
 delete []s;

 tama=0;

 for (int j=0;j<V.N_Antecedente(); j++){
     longitud = V.SizeDomain(j);
     if (mr[0][j]<mr[0][nr[0]-1]){
       for (int s=0; s<longitud; s++)
	 regla2[tama+s]='1';
     }
   tama+=longitud;
 }



  string patr3 = Union(regla1,regla2);
  aux =R1.Get_Code_Rule(r1);

  for (int j=0; j<nb[0]; j++)
    aux->PutValueBinary(0,j,patr3[j]);

}

string Union(const string &cad1, const string &cad2){
  string aux = cad1;

  for (int s=0; s<cad1.length(); s++)
    if (cad1[s]=='0' and cad2[s]=='0')
      aux[s]='0';
    else
      aux[s]='1';

  return aux;
}



// Hace la interseccion de las dos reglas y la operacion la almacena en la regla r1
void IntersecionReglasDiscriminantes(const VectorVar &V, ruleset &R1, int r1, ruleset &R2, int r2){
  genetcode r_1, r_2, *aux;
  aux=R1.Get_Code_Rule(r1);
  r_1=*aux;
  aux=R2.Get_Code_Rule(r2);
  r_2=*aux;

 int b,*nb,n,*nn,**nnn,r,*nr;
 char **mb;
 double **mr;
 string regla1, regla2, regla3;

 r_1.GetBinary(b,nb,mb);
 r_1.GetReal(r,nr,mr);

 char *s= new char[nb[0]+1];
 for (int j=0; j<nb[0]; j++)
   s[j]=mb[0][j];
 s[nb[0]]='\0';
 regla1=s;
 delete []s;

 int longitud, tama=0;

 for (int j=0;j<V.N_Antecedente(); j++){
     longitud = V.SizeDomain(j);
     if (mr[0][j]<mr[0][nr[0]-1]){
       for (int s=0; s<longitud; s++)
	 regla1[tama+s]='1';
     }
   tama+=longitud;
 }


 r_2.GetBinary(b,nb,mb);
 r_2.GetReal(r,nr,mr);

 s= new char[nb[0]+1];
 for (int j=0; j<nb[0]; j++)
   s[j]=mb[0][j];
 s[nb[0]]='\0';
 regla2=s;
 delete []s;

 tama=0;

 for (int j=0;j<V.N_Antecedente(); j++){
     longitud = V.SizeDomain(j);
     if (mr[0][j]<mr[0][nr[0]-1]){
       for (int s=0; s<longitud; s++)
	 regla2[tama+s]='1';
     }
   tama+=longitud;
 }



  string patr3 = Interseccion(regla1,regla2);
  aux =R1.Get_Code_Rule(r1);

  for (int j=0; j<nb[0]; j++)
    aux->PutValueBinary(0,j,patr3[j]);

  aux->GetReal(r,nr,mr);

  for (int j=0;j<V.N_Antecedente(); j++){
     if (mr[0][j]<mr[0][nr[0]-1]){
       mr[0][j]=mr[0][nr[0]-1];
     }
 }


  //cout << patr3 << endl << endl;

}


string Interseccion(const string &cad1, const string &cad2){
  string aux = cad1;

  for (int s=0; s<cad1.length(); s++)
    if (cad1[s]=='1' and cad2[s]=='1')
      aux[s]='1';
    else
      aux[s]='0';

  return aux;
}

int Grados_Libertad(string cad){
  int cont=0;
  for (int i=0; i<cad.length(); i++)
    if (cad[i]=='?')
      cont++;

  return cont;
}


void OrdenarVector(int *v, int tama){
  for (int i=0; i<tama; i++){
    for (int j=tama-1;j>i;j--)
      if (v[j]>v[i]){
	int aux=v[j];
	v[j]=v[i];
	v[i]=aux;
      }
  }
}


bool AnalizarConjuntosDeReglas(const VectorVar &V, ruleset &R, ruleset &R_Descrip, ruleset &R_mix, int nuevas){
  bool cambio=false;
  int *eliminables= new int[R.N_rule()];
  int n_eliminables=0;
  genetcode *aux_mix, *aux_R, *aux_R_des;
  int iguales, iguales2;
  bool *eliminar= new bool[R.N_rule()], eliminacion_retrasada;

  for (int i=0; i<R.N_rule(); i++){
    eliminar[i]=false;
  }

  for (int i=0; i<nuevas; i++){
    if (!eliminar[i]){
      eliminacion_retrasada=false;
      string res = Patron2Reglas(R_Descrip,R,V,i);
      string res3 = Patron1Regla(R_Descrip,V,i);
      string res4 = Patron1Regla(R,V,i);
      aux_mix =R_mix.Get_Code_Rule(i);
      int len = aux_mix->TamaBinary(0);
      char *cad = new char[len+1];
      for (int j=0; j<len; j++)
	cad[j]=aux_mix->GetValueBinary(0,j);
      cad[len]='\0';
      string patr = cad;
      delete [] cad;

      aux_R =R.Get_Code_Rule(i);
      len = aux_R->TamaBinary(0);
      cad = new char[len+1];
      for (int j=0; j<len; j++)
	cad[j]=aux_R->GetValueBinary(0,j);
      cad[len]='\0';
      string regla_R = cad;
      delete [] cad;

      aux_R_des =R_Descrip.Get_Code_Rule(i);
      len = aux_R_des->TamaBinary(0);
      cad = new char[len+1];
      for (int j=0; j<len; j++)
	cad[j]=aux_R_des->GetValueBinary(0,j);
      cad[len]='\0';
      string regla_R_des = cad;
      delete [] cad;

      if (res!=patr){
	cout << "No coinciden los patrones para la regla " << i << endl;
	cout << res << endl;
	cout << patr << endl;
	Pausa();
      }

      if (regla_R!=patr and Comparar2patrones(patr,regla_R)!=1){
	cout << "La regla discriminante "<< i << " no coincide con su patron\n";
	cout << regla_R << endl;
	cout << patr << endl;
	Pausa();
      }

      if (regla_R_des!=patr and Comparar2patrones(patr,regla_R_des)!=1){
	cout << "La regla descriptiva "<< i << " no coincide con su patron\n";
	cout << regla_R_des << endl;
	cout << patr << endl;
	Pausa();
     }

      //cout << "Regla " << i <<": " <<endl;
      //Comparar2patrones(res,patr);
      //cout << endl;
      for (int j=i+1; j<R.N_rule() and !eliminar[i]; j++){
	if (!eliminar[j]){
	  if (R.Class_of_rule(i)==R.Class_of_rule(j)){
	    string res2 = Patron2Reglas(R_Descrip,R,V,j);
	    //cout << "Reglas " << i <<" y " << j << " de la clase " << R.Class_of_rule(i) <<endl;
	    //Pinta2Reglas(R_Descrip,R,V,i);
	    //cout << endl;
	    //Pinta2Reglas(R_Descrip,R,V,j);
	    //cout << endl;

	    iguales = Comparar2patrones(res2,patr);
	    iguales2 = Comparar2patrones(res3,res2);
	    int iguales3 = Comparar2patrones(res4,res2);

	    if (iguales2!=2){
	      //cout << "\t\t\t\tCompatible la regla descriptiva " << i << " con el patron de la regla " << j << endl;
	    }

	    if (iguales3!=2){
	      //cout << "\t\t\t\tCompatible la regla discriminante " << i << " con el patron de la regla " << j << endl;
	    }


	    // El caso raro de campatibildades entre patrones
	    if (iguales>1 and iguales2!=2 and iguales3!=2){
	      cout << "Caso especial de incompatibilidades entre patrones y si entre reglas\n";
	      cout << res3 << endl;
	      cout << res4 << endl;
	      cout << patr << endl << endl;
	      cout << Patron1Regla(R_Descrip,V,j) << endl;
	      cout << Patron1Regla(R,V,j) << endl;
	      cout << res2 << endl << endl;
	      Pausa();
	    }


	    //cout << endl;

	    // Los dos patrones son iguales
	    if (iguales==0){
	      cambio=true;
	      double n_mas[2], n_menos[2];
	      eliminar[i]=true;
	      //cout << "\t\t\t\tCompatible la regla " << i << " con la regla " << j << endl;
	      //cout << "\t\t\t\t\tReglas " << i <<" y " << j << " de la clase " << R.Class_of_rule(i) <<endl;
	      //Pinta2Reglas(R_Descrip,R,V,i);
	      //cout << endl;
	      //Pinta2Reglas(R_Descrip,R,V,j);
	      //Pausa();
	    }
	    // Los dos patrones son compatibles
	    else if (iguales==1){
	      cambio=true;
	      cout << "\t\t\t\t\tCombinando Reglas " << i <<" y " << j << " de la clase " << R.Class_of_rule(i) <<endl;
		//cout << patr << endl;
		//cout << res2 << endl;
		//Comparar2patronesPinta(patr,res2);
		//cout << endl;
		//cout << res3 << endl;
		//cout << res4 << endl;
		//Comparar2patronesPinta(res3,res4);
		//Pausa();
	      double n_mas[4], n_menos[4];
	      n_mas[0] = R_Descrip.Get_N_Mas(i);
	      n_menos[0] = R_Descrip.Get_N_Menos(i);

	      n_mas[1] = R.Get_N_Mas(i);
	      n_menos[1] = R.Get_N_Menos(i);

	      n_mas[2] = R_Descrip.Get_N_Mas(j);
	      n_menos[2] = R_Descrip.Get_N_Menos(j);

	      n_mas[3] = R.Get_N_Mas(j);
	      n_menos[3] = R.Get_N_Menos(j);

	      //cout << "Regla "<< i << ": ";
	      //for (int t=0; t<2; t++)
		//cout << "("<< n_mas[t] << "/ "<< n_menos[t] << " = " << (n_mas[t]+1)/(n_mas[t]+n_menos[t]+1) << ") ";
	      //cout << endl;

	      //cout << "Regla "<< j << ": ";
	      //for (int t=2; t<4; t++)
		//cout << "("<< n_mas[t] << "/ "<< n_menos[t] << " = " << (n_mas[t]+1)/(n_mas[t]+n_menos[t]+1) << ") ";
	      //cout << endl;
	      //cout << "Patron regla : " << i << endl;
	      //Pinta2Reglas(R_Descrip,R,V,i);
	      string cad1 = Patron2Reglas(R_Descrip,R,V,i);
	      //cout << cad1 << endl;
	      //cout << "Patron regla : " << j << endl;
	      //Pinta2Reglas(R_Descrip,R,V,j);
	      string cad3 = Patron2Reglas(R_Descrip,R,V,j);
	      //cout << cad2 << endl;
	      //cout << "\t\t\t\tCombinadas las reglas " << i << " y " << j << endl;
	      //cout << "\t\t\t\t\tGrados de libertad antes   : " << Grados_Libertad(Patron2Reglas(R_Descrip,R,V,j)) << endl;

	      IntersecionReglasDiscriminantes(V,R,j,R,i);

	      UnionReglasDescriptivas(V,R_Descrip,j,R_Descrip,i);

	      string cad2 = Patron2Reglas(R_Descrip,R,V,j);


	      //cout << cad1 << endl;
	      //cout << cad3 << endl;
	      //cout << endl;
	      //cout << cad2 << endl;
	      string cad4 = Patron2Patrones(cad1,cad3);
	      cout <<  cad4 << endl;
	      int iguales4 = Comparar2patrones(cad2,cad4);
	      if (iguales4!=0){
		cout << cad2 << endl;
		cout <<  cad4 << endl;
		Comparar2patronesPinta(cad2,cad4);
		Pausa();
	      }

	      aux_mix =R_mix.Get_Code_Rule(j);
	      len = aux_mix->TamaBinary(0);
	      for (int t=0; t<len; t++)
		aux_mix->PutValueBinary(0,t,cad4[t]);

	      aux_R =R.Get_Code_Rule(j);
	      len = aux_R->TamaBinary(0);
	      for (int t=0; t<len; t++)
		if (cad4[t]=='?')
		  aux_R->PutValueBinary(0,t,'1');
		else
		  aux_R->PutValueBinary(0,t,cad4[t]);

	      aux_R_des =R_Descrip.Get_Code_Rule(j);
	      len = aux_R_des->TamaBinary(0);
	      for (int t=0; t<len; t++)
		if (cad4[t]=='?')
		  aux_R_des->PutValueBinary(0,t,'0');
		else
		  aux_R_des->PutValueBinary(0,t,cad4[t]);

	      if (Patron2Reglas(R_Descrip,R,V,j)!=cad4){
		cout << "Error al componer la descrip. y la predic. en la regla " << j << endl;
		Pausa();
	      }

	      //cout << "\t\t\t\t\tGrados de libertad despues : " << Grados_Libertad(cad4) << endl;

	      //Pausa();


	      eliminar[i]=true;
	      R_Descrip.Modificada(j,true);
	      R.Modificada(j,true);

	      if (Grados_Libertad(res2)>=Grados_Libertad(patr) ){
		R_Descrip.Add_N_Mas(j,n_mas[0]);
		R_Descrip.Add_N_Menos(j,n_menos[0]);

		R.Add_N_Mas(j,n_mas[1]);
		R.Add_N_Menos(j,n_menos[1]);
	      }

	      //R_Descrip.Add_N_Mas(j,n_mas[0]);
	      //R_Descrip.Add_N_Menos(j,n_menos[0]);

	      //R.Add_N_Mas(j,n_mas[1]);
	      //R.Add_N_Menos(j,n_menos[1]);

	      //n_mas[2] = R_Descrip.Get_N_Mas(j);
	      //n_menos[2] = R_Descrip.Get_N_Menos(j);

	      //n_mas[3] = R.Get_N_Mas(j);
	      //n_menos[3] = R.Get_N_Menos(j);

	      //cout << " Peso resultante: "<< "("<< n_mas[2] << "/ "<< n_menos[2] << " = " << (n_mas[2]+1)/(n_mas[2]+n_menos[2]+1) << ") ";
	      //cout << "                : "<< "("<< n_mas[3] << "/ "<< n_menos[3] << " = " << (n_mas[3]+1)/(n_mas[3]+n_menos[3]+1) << ")\n\n ";
	      //cout << "Nuevo Patron regla " << j << endl;
	      //cout << Patron2Reglas(R_Descrip,R,V,j) << endl;
	      //Pinta2Reglas(R_Descrip,R,V,j);
	      //PintaDiferenciaPatrones(cad2,Patron2Reglas(R_Descrip,R,V,j));

	      //char ch;
	      //cin >> ch;
	    }
	    else if (DistanciaHamming(res2,patr)==1){
	      //cout << res2 << endl;
	      //cout << patr << endl;
	      //PintaDiferenciaPatrones(res2,patr);
	      /*int pos = DiferenciaDistanciaHamming(res2,patr,0);
	      if (pos!=-1){
		res2[pos]='?';
		genetcode *ind = R_Descrip.Get_Code_Rule(j);
		ind->PutValueBinary(0,j,'0');
		ind = R.Get_Code_Rule(j);
		ind->PutValueBinary(0,j,'1');
		ind = R_mix.Get_Code_Rule(j);
		ind->PutValueBinary(0,j,'?');
		R_Descrip.Add_N_Mas(j,R_Descrip.Get_N_Mas(i));
		R_Descrip.Add_N_Menos(j,R_Descrip.Get_N_Menos(i));
		R.Add_N_Mas(j,R_Descrip.Get_N_Mas(i));
		R.Add_N_Menos(j,R_Descrip.Get_N_Menos(i));
		R_Descrip.Recalculate_Weight(j);
		R.Recalculate_Weight(j);
		eliminacion_retrasada=true;
	      }*/
	      //Pausa();
	    }
	  }
	  else if (R.Class_of_rule(i)!=R.Class_of_rule(j) and i<nuevas){
	    string res2 = Patron2Reglas(R_Descrip,R,V,j);
	    iguales = Comparar2patrones(res2,patr);

	    if (iguales==0){
	      cout << "las reglas " << i << " y " << j << " son contradictorias\n";
	      R.Add_N_Menos(i,R.Get_N_Mas(i));
	      R.Recalculate_Weight(i);
	      R_Descrip.Add_N_Menos(i,R_Descrip.Get_N_Mas(i));
	      R_Descrip.Recalculate_Weight(i);
	    }
	    else if (iguales==1){
	      cout << "las reglas " << i << " y " << j << " entran en conflicto\n";
	    }
	    else{
	      //no entra en principio en conflicto con ninguna otra regla
	    }

	  }
	}
      }
      if (eliminacion_retrasada)
	eliminar[i]=true;
    }
  }

  for (int i=0; i<R.N_rule()-1; i++){
    if (eliminar[i]){
      eliminables[n_eliminables]=i;
      n_eliminables++;
    }
  }


  //cout << "Pulse una tecla para continuar\n";
  //char ch;
  //cin >> ch;

  Filtrar_Reglas_New(R,R_Descrip,R_mix, eliminables, n_eliminables, false);


  delete [] eliminables;
  delete [] eliminar;

  return cambio;
}

bool Terminado_Aprendizaje(multipoblacion &G, example_set &E, bool &ultimo_modo){
  if ( (G.ValorFitness(0)>-1) or
        (  G.ValorFitness(0)<=-1 and G.ValorFitness(0,1)>1 and G.ValorFitness(0,1)>= (-G.ValorFitness(0))+2 ) or
        (G.ValorFitness(0)==-E.N_Examples() and G.ValorFitness(0,1)>1) or
       ((G.ValorFitness(0)==-E.N_Examples() ) and G.ValorFitness(0,1)==1 and G.ValorFitness(0,2)>0) or
       ( G.ValorFitness(0)> -E.N_Examples() and G.ValorFitness(0)<=-1) or
       ( (G.ValorFitness(0)==-E.N_Examples()) and G.ValorFitness(0,1)==1 and G.ValorFitness(0,2)==0 and G.ValorFitness(0,3)>0) or
       ( (G.ValorFitness(0)==-E.N_Examples()) and G.ValorFitness(0,1)==1 and G.ValorFitness(0,2)==0 and G.ValorFitness(0,3)==0 and G.ValorFitness(0,4)>0 )){
     ultimo_modo=false;
     return false;
   }
   else if (!ultimo_modo){
          cout <<"\t================================================\n";
          cout <<"\t=                                              =\n";
          cout <<"\t=                 MODO FINAL                   =\n";
          cout <<"\t=                                              =\n";
          cout <<"\t================================================\n";
          ultimo_modo=true;
          return false;
        }
        else
          return true;

}


bool Cambios(int n_act, double **m_act, int n_old, double **m_old, int cols)
{
    bool iguales=true;

    if (n_act!=n_old)
        return true;
    else
    {
        int i=0,j;
        while (i<n_act && iguales)
        {
            j=0;
            while (j<cols && m_act[i][j]==m_old[i][j])
                j++;
            if (j==cols)
                i++;
            else
                iguales=false;
        }
        return !iguales;
    }

}



int ** Encontrar_Individuos_dados_ejemplos(const VectorVar &V, example_set &E, int nclases, int n_individuos, double ***Tabla, int *n_ejemplos_por_clase){
  int a;

  int **matriz = new int*[n_individuos];
  for (int i=0; i< n_individuos; i++)
    matriz[i] =new int[V.N_Antecedente()];

  for (int i=0; i<n_individuos; i++){
    // Elegir un ejemplo de la clase i%n_clases
    if (n_ejemplos_por_clase[i%nclases]==0){
      for (int j=0; j<V.N_Antecedente(); j++)
         matriz[i][j]=-1;
    }
    else {
       a=Select_Random_Individuo(E.N_Examples(),-1,0);
       while (E.Is_Covered(a) || T_Adaptacion_Variable_Etiqueta(V,Tabla,a,V.Consecuente(),i%nclases)==0){
         a++;
         a=a%E.N_Examples();
       }
       for (int j=0; j<V.N_Antecedente(); j++){
         matriz[i][j]=BetterEtiqueta(V,Tabla,j,a);
       }
    }
  }

  return matriz;
}




int ** Encontrar_Individuos_dados_ejemplos_para_clase(const VectorVar &V, example_set &E, int nclases, int n_individuos, double ***Tabla, int *n_ejemplos_por_clase, int clase){
  int a;

  int **matriz = new int*[3];
  for (int i=0; i< 3; i++)
    matriz[i] =new int[V.N_Antecedente()];

  for (int i=0; i<3; i++){
    // Elegir un ejemplo de la clase "clase"
    if (n_ejemplos_por_clase[clase]==0){
      for (int j=0; j<V.N_Antecedente(); j++)
         matriz[i][j]=-1;
    }
    else {
       a=Select_Random_Individuo(E.N_Examples(),-1,0);
       while (E.Is_Covered(a) || T_Adaptacion_Variable_Etiqueta(V,Tabla,a,V.Consecuente(),clase)==0){
         a++;
         a=a%E.N_Examples();
       }
       for (int j=0; j<V.N_Antecedente(); j++){
         matriz[i][j]=BetterEtiqueta(V,Tabla,j,a);
       }
    }
  }

  return matriz;
}



void Obtener_Aciertos_y_Errores_Reglas_Aprendidas(example_set &E, int *basica_para_regla, int *errores_para_regla, int nrule){

  for (int i=0; i<nrule; i++){
    basica_para_regla[i]=0;
    errores_para_regla[i]=0;
  }

  for (int j=0; j<E.N_Examples(); j++) {
     if (E.Rule_Aislada(j)==1 and E.Rule_Used(j)>=0){
        basica_para_regla[E.Rule_Used(j)]++;
     }
     else if (E.Rule_Aislada(j)==-1 and E.Rule_Used(j)>=0){
             errores_para_regla[E.Rule_Used(j)]++;
          }
  }
}




void ReajustarPesoaReglas(multipoblacion &G, const VectorVar &V, example_set &E, ruleset &R, VectorVar &F,
       double ***Tabla, double ***Tabla_Nuevas){

  genetcode *aux;
  double n_mas, n_menos;
  for (int i=0; i<R.N_rule(); i++){
    //if (R.Get_Modificada(i)){
    if (true){
      cout << "Peso antes de Ajustar: (" << R.Get_N_Mas(i) << "/ "<< R.Get_N_Menos(i) << ") = " << R.Get_Weight(i) << endl;
      R.Anular_N_Mas_y_N_Menos(i);
      aux=R.Get_Code_Rule(i);
      G.PutCodigo(G.N_individuos()-1,*aux);
      double peso = R.Get_Weight(i);
      cout << "Rule " << i << " pasa de (" << R.Get_N_Mas(i) << "/ "<< R.Get_N_Menos(i) << ") = " << peso ;
      peso = Weight_Individual(G,V,E,F,Tabla, Tabla_Nuevas, G.N_individuos()-1, n_mas, n_menos);
      R.Add_Weight(i,n_mas, n_menos);
      cout << " a (" << R.Get_N_Mas(i) << "/ "<< R.Get_N_Menos(i) << ") = " << peso << endl;
      R.Modificada(i,false);
    }
  }

  /*char ch;
  cin >> ch;

  cout << " Proceso de adaptacion de las reglas existentes a los nuevos ejemplos\n";
  for (int i=0; i<R.N_rule(); i++){
    aux=R.Get_Code_Rule(i);
    for (int j=0; j<G.N_individuos; j++)
      G.PutCodigo(j,*aux);
  }*/
}


void ActualizarPeso(multipoblacion &G, const VectorVar &V, example_set &E, ruleset &R, VectorVar &F,
       double ***Tabla, double ***Tabla_Nuevas){

  genetcode *aux;
  double n_mas, n_menos;
  for (int i=0; i<R.N_rule(); i++){
    if (R.Get_Modificada(i)){
    //if (true){
      //cout << "Peso antes de Ajustar: (" << R.Get_N_Mas(i) << "/ "<< R.Get_N_Menos(i) << ") = " << R.Get_Weight(i) << endl;
      aux=R.Get_Code_Rule(i);
      G.PutCodigo(G.N_individuos()-1,*aux);
      double peso = R.Get_Weight(i);
      //cout << "Rule " << i << " pasa de (" << R.Get_N_Mas(i) << "/ "<< R.Get_N_Menos(i) << ") = " << peso ;
      peso = Weight_Individual(G,V,E,F,Tabla, Tabla_Nuevas, G.N_individuos()-1, n_mas, n_menos);
      R.Add_Weight(i,n_mas, n_menos);
      //cout << " a (" << R.Get_N_Mas(i) << "/ "<< R.Get_N_Menos(i) << ") = " << peso << endl;
      R.Modificada(i,false);
    }
  }
}




int GA(multipoblacion &G, const VectorVar &V, example_set &E, ruleset &R, ruleset &R_dis, ruleset &R_mix, double **I, double *agregado,
       double *peso_agregado, int Total_ejemplos, bool discriminante, bool primeraVez, bool reducir_reglas,
       int *ponderacion_por_clase, double ***Tabla, double ***Tabla_Nuevas, VectorVar &F, int *basica_para_regla, int *errores_para_regla,
       int *reglas_a_eliminar, int &n_reglas_a_eliminar, bool &Added_Rule, bool &Specific_Rule, int tam_rule_not_removed)
{

    Added_Rule = false;

    double peso;
    int NFitness=G.N_Valoracion();
    double *valor = new double [NFitness];
    int n_clases=V.SizeDomain(V.Consecuente());
    int *n_ejemplos_por_clase = new int[n_clases];
    ofstream f;
    int n=G.N_individuos();
    int nbetter; // Numero de los mejores individuos
    int old_nbetter; // Numero de los mejores en t-1
    int *better = new int[n]; // Lista de los mejores individuos
    double **fbetter = new double*[n]; // Valoracion de los mejores individuos
    for (int i=0; i<n; i++)
        fbetter[i]=new double[NFitness];
    double **fbetter_old = new double*[n]; // Valoracion de los mejores individuos en t-1
    for (int i=0; i<n; i++)
        fbetter_old[i]=new double[NFitness];

    int sin_cambios = 0;
    int z;
    double aciertos;


  E.Examples_per_Class(V.Consecuente(),0,n_clases, n_ejemplos_por_clase);
  int **individuos_por_clase = new int*[n_clases];
  for (int i=0; i<n_clases; i++){
    individuos_por_clase[i] = new int[n+1];
    individuos_por_clase[i][0]=0;
  }


  int **sujetos, *lista_tamano;

  if (R.N_rule()==0 or primeraVez){

    sujetos = Encontrar_Individuos_dados_ejemplos(V,E,n_clases,G.N_individuos(),Tabla, n_ejemplos_por_clase);
    lista_tamano = new int[V.N_Antecedente()];
    for (int i=0; i<V.N_Antecedente(); i++)
      lista_tamano[i]=V.SizeDomain(i);

    G.PoblacionInicial(I, n_clases, E.N_Examples(), sujetos, V.N_Antecedente(), lista_tamano); // Revisar

    for (int i=0; i<G.N_individuos(); i++)
      delete [] sujetos[i];
    delete [] sujetos;
    delete [] lista_tamano;
  }
  else {
    /*sujetos = Encontrar_Individuos_dados_ejemplos_para_clase(V,E,n_clases,G.N_individuos(),Tabla, n_ejemplos_por_clase, G.ClaseIndividuo(0));
    lista_tamano = new int[V.N_Antecedente()];
    for (int i=0; i<V.N_Antecedente(); i++)
      lista_tamano[i]=V.SizeDomain(i);

    G.PoblacionInicial2_Estacionario(I, n_clases, E.N_Examples(), sujetos, V.N_Antecedente(), lista_tamano, G.ClaseIndividuo(0));*/


    for (int i=0; i<G.N_individuos(); i++)
      G.PutModificado(i);

    /*for (int i=0; i<3; i++)
      delete [] sujetos[i];
    delete [] sujetos;
    delete [] lista_tamano;*/
    }

    //G.PoblacionInicialAleatoria_4L(I, n_clases, E.N_Examples());



    fitness(G,V,E,R, discriminante, reducir_reglas ,ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F, basica_para_regla, errores_para_regla, tam_rule_not_removed);

    genetcode code;
    G.Codigo(G.N_individuos()-1,code);
    int len = code.TamaBinary(0);
    char *cad = new char[len+1];
    for (int i=0; i<len; i++)
      cad[i]='?';
    cad[len]='\0';
    string patron_general = cad;
    delete [] cad;
    string vacia="";

    for (int ik=G.N_individuos()-2; ik<G.N_individuos(); ik++)
      Obtener_Regla_Descriptiva(G,V,E,R,ik,discriminante,reducir_reglas,ponderacion_por_clase,Tabla, Tabla_Nuevas, F, basica_para_regla, errores_para_regla,
					    reglas_a_eliminar, n_reglas_a_eliminar, false, patron_general, true, vacia, tam_rule_not_removed);

    //G.Sort_4L(E, n_clases, n_ejemplos_por_clase, individuos_por_clase, agregado, Total_ejemplos);
    //G.Sort(n_clases, n_ejemplos_por_clase);
    G.Sort();
    G.Better(n_clases, better, old_nbetter, individuos_por_clase);
    for (z=0; z<old_nbetter; z++)
        for (int l=0; l<NFitness; l++)
            fbetter_old[z][l]=G.ValorFitness(better[z],l);

    int t=0;
    //f.open("slave.log",ios::out);
    //f << "ITERACION " << t << endl;
    //f.close();
    //G.PintaInFile();

    int max_iter=500;


    if (G.ValorFitness(0)!=-E.N_Examples() )
       max_iter=500;




    while (/*t<=1000*/ sin_cambios<=max_iter)
    {
        //G.Seleccion();
        //G.Cruce2puntos_Estacionario();
        if ( (R.N_rule()>1 and G.ValorFitness(0)<=0 and Probabilidad(0.75)) or (R.N_rule()>1 and Probabilidad(0.01))){
          int pos=Select_Random_Individuo(R.N_rule(),-1,-1);
          //cout << ">>>>>>>>>>>>Regla Insertada: " << pos << endl;
          genetcode *aux;
          /*for (int k=G.N_individuos()-2; k<G.N_individuos();k++){

	    if (k==G.N_individuos()-2)
	      aux=R.Get_Code_Rule(pos);
	      //Obtener_Regla_Descriptiva(G,V,E,R,k,discriminante,reducir_reglas,ponderacion_por_clase,Tabla, Tabla_Nuevas, F, basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar, false, patr);
	    else
	      aux=R_dis.Get_Code_Rule(pos);
	      //Obtener_Regla_Discriminante(G,V,E,R,k,discriminante,reducir_reglas,ponderacion_por_clase,Tabla, Tabla_Nuevas, F, basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar, false, patr);

            G.PutCodigo(k,*aux);
            pos=Select_Random_Individuo(R.N_rule(),pos,-1);
            //cout << ">>>>>>>>>>>>Regla Insertada: " << pos << endl;
	    //G.PutModificado(k);
          }*/
          G.MutacionUniforme_Estacionario_Especial();
        }
        else {
          G.CruceBasedLogical_Estacionario((1.0*sin_cambios)/max_iter);
          G.MutacionUniforme_Estacionario();
        }
        fitness(G,V,E, R, discriminante, reducir_reglas ,ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F, basica_para_regla, errores_para_regla, tam_rule_not_removed);
	aciertos=G.Sort_4L(E, n_clases, n_ejemplos_por_clase, individuos_por_clase, agregado, peso_agregado, Total_ejemplos);
        G.Better(n_clases, better, nbetter, individuos_por_clase);
        for (z=0; z<nbetter; z++)
        {
            for (int l=0; l<NFitness; l++)
                fbetter[z][l]=G.ValorFitness(better[z],l);
        }


        if (sin_cambios==0)
        {
            //cout << t <<": " << endl;
            //G.PintaFitness_Estacionario(0);
            for (z=0; z<nbetter; z++)
            {
                //cout <<"\t";
                //G.PintaFitness_Estacionario(better[z]);
		//G.PintaIndividuo(0);
		/*genetcode code;
		G.Codigo(0,code);
		R.Add(code,peso);
		cout << R.SPrint(R.N_rule()-1) << endl;
		R.Remove();*/
                //cout << endl;
            }
            //cout <<"Aciertos: "<< aciertos << "  ("<< 1.0*aciertos/Total_ejemplos <<")" <<endl;
            //char ch;
            //cin >> ch;
        // Controlar numero de individuos de cada clase en la poblacion final
        /*int *ind_por_clase = new int[n_clases];
        double *mejor_de_clase = new double[n_clases];

        for (int l=0; l<n_clases; l++){
          ind_por_clase[l]=0;
          mejor_de_clase[l]=-E.N_Examples();
        }

        for (int l=0; l<G.N_individuos();l++){
          if (mejor_de_clase[G.ClaseIndividuo(l)]<G.ValorFitness(l))
            mejor_de_clase[G.ClaseIndividuo(l)]=G.ValorFitness(l);
          ind_por_clase[G.ClaseIndividuo(l)]++;
        }

        cout << "\t\t\t\t " << ind_por_clase[0] << "("<<mejor_de_clase[0]<<")";
        for (int l=1; l<n_clases; l++)
          cout << ", " << ind_por_clase[l]<< "("<<mejor_de_clase[l]<<")";
        cout << endl;

        delete [] ind_por_clase;
        delete [] mejor_de_clase;*/
        }
        t++;

        //f.open("slave.log",ios::app);
        //f << "ITERACION " << t << endl;
        //f.close();
        //G.PintaInFile();

        if (!Cambios(nbetter, fbetter, old_nbetter, fbetter_old, NFitness)){
            sin_cambios++;
	    if (sin_cambios==2){
	      int t;
	      //t= G.N_individuos()-1;
	      t=1;
	      genetcode code;
	      G.Codigo(0,code);
	      G.PutCodigo(t,code);


	      int len = code.TamaBinary(0);
	      char *cad = new char[len+1];
	      for (int i=0; i<len; i++)
		cad[i]='?';
	      cad[len]='\0';
	      string patr = cad;
	      delete [] cad;
	      string cadena_vacia="";

	      if (G.ValorFitness(t)!=-E.N_Examples() ){
		  Obtener_Regla_Discriminante(G,V,E,R,t,discriminante, reducir_reglas,ponderacion_por_clase,Tabla, Tabla_Nuevas, F, basica_para_regla, errores_para_regla,
					    reglas_a_eliminar, n_reglas_a_eliminar, false, patr, true, cadena_vacia, tam_rule_not_removed);
	      }
	      else {
		if (G.ValorFitness(t,1)==1){
		  Obtener_Regla_Descriptiva(G,V,E,R,t,discriminante,reducir_reglas,ponderacion_por_clase,Tabla, Tabla_Nuevas, F, basica_para_regla, errores_para_regla,
					    reglas_a_eliminar, n_reglas_a_eliminar, false, patr, true, cadena_vacia, tam_rule_not_removed);
		}
		else {
		  Obtener_Regla_Discriminante(G,V,E,R,t,discriminante,reducir_reglas,ponderacion_por_clase,Tabla, Tabla_Nuevas, F, basica_para_regla, errores_para_regla,
					    reglas_a_eliminar, n_reglas_a_eliminar, false, patr, true, cadena_vacia, tam_rule_not_removed);
		}
	      }
	    }
	}
        else
        {
            sin_cambios=0;

	    if (G.ValorFitness(0)!=-E.N_Examples() ){
		max_iter=500;
	    }

	    for (z=0; z<nbetter; z++)
		  for (int l=0; l<NFitness; l++)
		      fbetter_old[z][l]=fbetter[z][l];
	      old_nbetter=nbetter;
        }

    }

    // Controlar numero de individuos de cada clase en la poblacion final
    /*int *ind_por_clase = new int[n_clases];
    double *mejor_de_clase = new double[n_clases];

    for (int l=0; l<n_clases; l++){
      ind_por_clase[l]=0;
      mejor_de_clase[l]=-E.N_Examples();
    }

    for (int l=0; l<G.N_individuos();l++){
          if (mejor_de_clase[G.ClaseIndividuo(l)]<G.ValorFitness(l))
            mejor_de_clase[G.ClaseIndividuo(l)]=G.ValorFitness(l);
      ind_por_clase[G.ClaseIndividuo(l)]++;
    }

    cout << "=====\t\t " << ind_por_clase[0] << "("<<mejor_de_clase[0]<<")";
    for (int l=1; l<n_clases; l++)
       cout << ", " << ind_por_clase[l]<< "("<<mejor_de_clase[l]<<")";
    cout << endl;
    //cin.get();

    delete [] ind_por_clase;
    delete [] mejor_de_clase;*/

    // Para sacar la mejor regla de cada clase

    //cout << "--------------------------------------------------"<<endl;
    G.Better(n_clases, better, nbetter, individuos_por_clase);

    delete [] n_ejemplos_por_clase;

    if (  G.ValorFitness(0)<=-1 and G.ValorFitness(0,1)>1 and G.ValorFitness(0,1)>= (-G.ValorFitness(0))+2 ){

      cout << "->->->->->-> El nuevo caso especial\n";
      //sleep(2);
    }

    //Bondad_Final(G,V,E,R,valor, NFitness,0, discriminante, reducir_reglas, ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F, basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar);
   //G.Valoracion(0,valor);


   //cout << "------------------------>>>>>>> Mejor Individiuo\n";
   //G.PintaFitness_Estacionario(0);
   //cout << "------------------------>>>>>>>*****************\n";


    if ((G.ValorFitness(0)>-1) or
        (  G.ValorFitness(0)<=-1 and G.ValorFitness(0,1)>1 and G.ValorFitness(0,1)>= (-G.ValorFitness(0))+2 ) or
        (  G.ValorFitness(0)==-E.N_Examples() and G.ValorFitness(0,1)>1) or
        ( (G.ValorFitness(0)==-E.N_Examples()) and G.ValorFitness(0,1)==1 and G.ValorFitness(0,2)>0) or
        ( (G.ValorFitness(0)==-E.N_Examples()) and G.ValorFitness(0,1)==1 and G.ValorFitness(0,2)==0 and G.ValorFitness(0,3)>0) or
        ( (G.ValorFitness(0)==-E.N_Examples()) and G.ValorFitness(0,1)==1 and G.ValorFitness(0,2)==0 and G.ValorFitness(0,3)==0 and G.ValorFitness(0,4)>0 )){
      genetcode code, code1, code2;
      double val_aux = G.ValorFitness(0);
      G.Codigo(0,code);

      int len = code.TamaBinary(0);
      char *cad = new char[len+1];
      for (int i=0; i<len; i++)
	cad[i]='?';
      cad[len]='\0';
      string patr = cad;
      delete [] cad;


      double valor2[10], nmas1, nmas2, nmenos1, nmenos2, esp[4], gen[4];


      //Aciertos(R,V,E, F);


      string s_esp,s_gen, cadena_vacia="";

      do {
	//R.Anular_N_Mas_y_N_Menos(r_obj);
	//R_dis.Anular_N_Mas_y_N_Menos(r_obj);

	Obtener_Regla_Discriminante(G,V,E,R,0,discriminante, reducir_reglas,ponderacion_por_clase,Tabla,
				    Tabla_Nuevas, F, basica_para_regla, errores_para_regla,
				    reglas_a_eliminar, n_reglas_a_eliminar, false, patr, true, cadena_vacia, tam_rule_not_removed);

        G.Codigo(0,code1);
	len = code1.TamaBinary(0);
	cad = new char[len+1];
	for (int i=0; i<len; i++)
	  cad[i]=code1.GetValueBinary(0,i);
	cad[len]='\0';
	s_gen = cad;
	delete [] cad;



	Bondad(G,V,E,R,valor,NFitness,0, discriminante, reducir_reglas ,ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
		 basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar, false, cadena_vacia, tam_rule_not_removed);

        peso = Weight_Individual(G,V,E,F,Tabla, Tabla_Nuevas, 0, nmas1, nmenos1);

	gen[0]=valor[0];
	gen[1]=valor[1];
	gen[2]=0;
	gen[3]=0;

	Obtener_Regla_Descriptiva(G,V,E,R,0,discriminante, reducir_reglas,ponderacion_por_clase,Tabla,
				    Tabla_Nuevas, F, basica_para_regla, errores_para_regla,
				    reglas_a_eliminar, n_reglas_a_eliminar, false, patr, true, cadena_vacia, tam_rule_not_removed);

	Bondad(G,V,E,R,valor,NFitness,0, discriminante, reducir_reglas ,ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
		 basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar, false, cadena_vacia, tam_rule_not_removed);

        peso = Weight_Individual(G,V,E,F,Tabla, Tabla_Nuevas, 0, nmas2, nmenos2);

	esp[0]=valor[0];
	esp[1]=valor[1];
	esp[2]=0;
	esp[3]=0;

	G.Codigo(0,code2);
	len = code2.TamaBinary(0);
	cad = new char[len+1];
	for (int i=0; i<len; i++)
	  cad[i]=code2.GetValueBinary(0,i);
	cad[len]='\0';
	s_esp = cad;
	delete [] cad;



	/*if ( !ComprobacionEspGen(s_esp,s_gen) ){
	  cout << "esp: " << s_esp << endl;
	  cout << "gen: " << s_gen << endl;
	  cout << "pat: ";
	  Comparar2patronesPinta(s_esp,s_gen);
	  cout << endl;
	}*/

	//Evento1_Generalizar(G, V, E, R, R_dis, R_mix, ponderacion_por_clase, F, Tabla, Tabla_Nuevas,r_obj,gen);
	//Evento1_Especificar(G, V, E, R, R_dis, R_mix, ponderacion_por_clase, F, Tabla, Tabla_Nuevas,r_obj,esp);
	//Pausa();
      } while (!vectoresIguales(esp,gen,4) or !ComprobacionEspGen(s_esp,s_gen));



      //cout << "Introduciendo reglas\n";

      R.Add(code2, nmas2, nmenos2, false);

      //cout << "-------------------------------------------------> Metida la regla\n";
      //cout << R.SPrint(R.N_rule()-1) << endl;


      R_dis.Add(code1, nmas1, nmenos1, false);

      patr = Patron2Reglas(R,R_dis,V,R.N_rule()-1);


      if ( !ComprobacionEspGen(Patron1Regla(R,V,R.N_rule()-1), Patron1Regla(R_dis,V,R_dis.N_rule()-1)) ){
	cout << "esp: " << endl << Patron1Regla(R,V,R.N_rule()-1) << endl;
	cout << Patron1Regla(R_dis,V,R_dis.N_rule()-1) << endl <<":gen" << endl;
	Comparar2patronesPinta(Patron1Regla(R,V,R.N_rule()-1), Patron1Regla(R_dis,V,R_dis.N_rule()-1));
	Pausa();
      }


      for (int i=0; i<len; i++)
	code1.PutValueBinary(0,i,patr[i]);
      R_mix.Add(code1,peso);



      /*code2=(*R_mix.Get_Code_Rule(R_mix.N_rule()-1));
      cout << "Patron de la regla obtenida: \n";
      for (int i=0; i<len; i++)
	cout << code2.GetValueBinary(0,i);
      cout << endl;

      Pausa();*/

      /*Obtener_Aciertos_y_Errores_Reglas_Aprendidas(E, basica_para_regla, errores_para_regla, R.N_rule());

      Bondad(G,V,E,R,valor,NFitness,0, discriminante, reducir_reglas ,ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
		 basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar, false, -1);
      valor[0]=esp[0];
      valor[1]=valor[1]-1;
      n_reglas_a_eliminar--;
      G.Valoracion(0,valor);

      G.PintaFitness_Estacionario(0);*/

      Added_Rule=true;

      if (  G.ValorFitness(0)==E.N_Examples() and G.ValorFitness(0,1)==1)
	Specific_Rule=true;
      else
	Specific_Rule=false;

      //cout << "------------------------>>>>>>> Regla Aniadida\n";
      //G.PintaFitness_Estacionario(0);
      //cout << "------------------------>>>>>>>>>>>>>>>>>>>>>>\n";
     /* cout << "-------------------------------------------------> Metida la regla\n";
      cout << "-------------------------------------------------> Descriptiva\n";
      cout << R.SPrint(R.N_rule()-1) << endl;
      cout << "\n-------------------------------------------------> Discriminante\n";
      cout << R_dis.SPrint(R.N_rule()-1) << endl;*/
      z=1;
      // Meter mas reglas
      //G.PutCodigo(0,code2);

      // Pongo la regla por defecto para la clase obtenida en la posición penultima de la poblacion
      for (int i=0; i<len; i++)
				code1.PutValueBinary(0,i,'1');

      int c = code1.GetValueInteger(0,0);

      G.PutCodigo(0,code2);
      code2.PutValueInteger(0,0, (c+1)%n_clases);

      //G.PutCodigo(G.N_individuos()-2,code1);


      //Pongo la regla por defecto para otra clase en la ultima posicion de la poblacion
      code1.PutValueInteger(0,0, (c+1)%n_clases);
      G.PutCodigo(G.N_individuos()-1,code1);
      c=(c+1)%n_clases;
      code1.PutValueInteger(0,0, (c+1)%n_clases);
      G.PutCodigo(G.N_individuos()-2,code1);


      while (z<nbetter)
      {
        G.Codigo(better[z],code);
        peso = Bondad(G,V,E,R,valor,NFitness,better[z],discriminante,reducir_reglas, ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
		      basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar, false, cadena_vacia, tam_rule_not_removed);
        R.Add(code,peso);
        z++;
      }
    }

    for (int i=0; i<n_clases; i++)
    {
        delete [] individuos_por_clase[i];
    }
    delete [] individuos_por_clase;
    delete [] better;
    for (int i=0; i<n; i++)
    {
        delete [] fbetter[i];
        delete [] fbetter_old[i];
    }
    delete [] fbetter;
    delete [] fbetter_old;
    delete [] valor;

    return t;
}



// Determina si se mantiene que esp es una especificacion de gen y viceversa
bool ComprobacionEspGen(string esp, string gen){
  bool correcto=true;

  for (int i=0; i<esp.length() and correcto; i++)
    if (esp[i]=='1' and gen[i]=='0')
      correcto=false;

  return correcto;
}



void fitness (multipoblacion &G, const VectorVar &V, const example_set &E, ruleset &R, bool discriminante,
              bool reducir_reglas, int *ponderacion_por_clase, int nrule, double ***Tabla, double ***Tabla_Nuevas, VectorVar &F,
	      int *basica_para_regla, int *errores_para_regla, int tam_rule_not_removed)
{

    string cadena_vacia="";
    int n=G.N_individuos();
    int n_fitnes_componentes=G.N_Valoracion();
    double *valor = new double [n_fitnes_componentes];
    int *reglas_a_eliminar = new int [R.N_rule()], n_reglas_a_eliminar;

		#pragma omp parallel for
    for (int i=0; i<n; i++)
    {
        if (G.Modificado(i))
        {
            Bondad(G,V,E,R,valor,n_fitnes_componentes,i, discriminante, reducir_reglas ,ponderacion_por_clase, nrule, Tabla, Tabla_Nuevas, F,
			 basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar, false, cadena_vacia, tam_rule_not_removed);
            G.Valoracion(i,valor);
        }
    }

    delete [] valor;
    delete [] reglas_a_eliminar;

}


bool Relevante( char *s, int tama){
 int i=0;
 while (i<tama and s[i]=='1')
   i++;

 return (i<tama);
}




bool Asignacion_Comprensible (char *s, int tama)
{
    int sec_unos=0, sec_ceros=0, n_ceros=0;
    bool act_unos = (s[0]=='1');

    if (act_unos)
        sec_unos=1;
    else
    {
        sec_ceros=1;
        n_ceros=1;
    }

    for (int i=1; i<tama; i++)
    {
        if (!act_unos && s[i]=='1')
        {
            act_unos=true;
            sec_unos++;
        }
        else if (!act_unos && s[i]=='0')
        {
            n_ceros++;
        }
        else if (act_unos && s[i]=='0')
        {
            sec_ceros++;
            n_ceros++;
            act_unos=false;
        }
    }

    return ( (sec_unos==1) );
}


int NumeroCeros (char *s, int tama)
{
    int acu=0;
    for (int i=0; i<tama; i++)
        if (s[i]=='0')
            acu++;
    return acu;
}



int NumeroUnos (char *s, int tama)
{
    int acu=0;
    for (int i=0; i<tama; i++)
        if (s[i]=='1')
            acu++;
    return acu;
}



void Obtener_Regla_Discriminante(multipoblacion &G, const VectorVar &V, example_set &E, ruleset &R, int elemento_i,
                               bool discriminante, bool reduccion_reglas, int* ponderacion_por_clase, double ***Tabla,
                               double ***Tabla_Nuevas, VectorVar &F, int *basica_para_regla,
			       int *errores_para_regla, int *reglas_a_eliminar, int &n_reglas_a_eliminar, bool debug, string patr, bool inGA,
			       string initial_rule, int tam_rule_not_removed)
{
    int r1;
    double *nr1=G.Subpoblacion_Real(0,elemento_i,r1);
    double umbral = nr1[r1-1], aux;
    int N_valor=G.N_Valoracion();
    double *valor = new double [N_valor], *valor2 = new double[N_valor];
    string cadena;
    bool mejora=false, mejora_global=false;
    int tama=0;
    int b1;
    char *nb1=G.Subpoblacion_Binaria(0,elemento_i,b1);
    char *s= new char[b1+1];
    for (int i=0; i<b1; i++)
      s[i]=nb1[i];
    s[b1]='\0';
    cadena = s;
    delete [] s;
    int longitud;
    for (int j=0; j<r1-1; j++)
    {
      longitud = V.Variable(j).N_etiquetas();
      if (nr1[j]<umbral){
	for (int s=0; s<longitud; s++){
	  cadena[tama+s]='1';
	  nb1[tama+s]='1';
	}
      }
      tama+=longitud;
    }

    //Aciertos(R,V,E,F);

    Bondad(G,V,E,R,valor,N_valor,elemento_i,discriminante, reduccion_reglas, ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
          basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar,!inGA, initial_rule, tam_rule_not_removed);


    G.Valoracion(elemento_i,valor);
    //cout << cadena << endl;
    //cout << "Valor antes de Discriminacion:   ";
    if (debug){
      cout << "Valor antes de Discriminacion:   " << endl;
      cout << cadena << endl;
      cout << "Valor de partida\n";
      G.PintaFitness_Estacionario(elemento_i);
    }

    /*for (int i=0; i<r1; i++)
      cout << "\t" << nr1[i];
    cout << endl;*/
    //cout << "Umbral: " << umbral << endl;

    /*if (umbral<=0){
      umbral==0.01;
      nr1[r1-1]=umbral;
    }*/

    tama=0;

    bool *se_puede_anular = new bool[r1-1];

    for (int j=0; j<r1-1; j++)
    {
      longitud = V.Variable(j).N_etiquetas();
      if (nr1[j]>=umbral)
	se_puede_anular[j]=true;
      else
	se_puede_anular[j]=false;

      for (int s=0; s<longitud and se_puede_anular[j]; s++)
	if (patr[tama+s]=='0')
	  se_puede_anular[j]=false;

      if (patr[tama]=='x')
	nr1[j]=umbral-0.001;

      tama+=longitud;
    }

    for(int j=0; j<patr.length(); j++)
      if (patr[j]=='0' or patr[j]=='1')
	nb1[j]=patr[j];


    tama=0;
    for (int j=0; j<r1-1; j++)
    {
      longitud = V.Variable(j).N_etiquetas();

        if (nr1[j]>=umbral and se_puede_anular[j])
        {
            aux = nr1[j];
	    mejora=false;
            nr1[j]=umbral-0.001;

	    Bondad(G,V,E,R,valor2,N_valor,elemento_i,discriminante, reduccion_reglas, ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
		  basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar,!inGA, initial_rule, tam_rule_not_removed);
	    //cout << "\tValor[0]= " << valor2[0] << endl;
		  if ( (valor2[0]>valor[0] /*and valor2[0]>0*/) or
		       (valor2[0]==valor[0] and valor2[1]>valor[1]) or
		       (valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]>valor[2]) /*or
		       (valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]==valor[2] and valor2[3]>=valor[3]) or
		       (valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]==valor[2] and valor2[3]==valor[3] and valor2[4]>=valor[4])*/
		    )
		  {
		  valor[0]=valor2[0];
		  valor[1]=valor2[1];
		  valor[2]=valor2[2];
		  valor[3]=valor2[3];
		  valor[4]=valor2[4];
		  mejora=true;
		  mejora_global=true;
		  //G.Valoracion(elemento_i,valor2);
		  //cout << "\t\t\t" << nb1 << endl;
		  //cout << "Anulada la variable " << V.Variable(j).Name() << "  Valor: " << valor2[0] << endl;
		  if (debug){
		    //cout << "\t\t\t" << nb1 << endl;
		    cout << "Anulada la variable " << V.Variable(j).Name() << "  Valor: " << valor2[0] << endl;
		    G.Valoracion(elemento_i,valor2);
		    G.PintaFitness_Estacionario(elemento_i);
		  }
	    }
	    else{
	      nr1[j]=aux;
	      //Bondad(G,V,E,R,valor2,N_valor,elemento_i,discriminante, true, ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
		  //basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar);
	      //cout << "No se puede anular la variable " << V.Variable(j).Name() << "  Valor: " << valor2[0] << endl;
	    }

        }


      tama+=longitud;
    }

    tama=0;
    for (int j=0; j<r1-1; j++)
    {
      longitud = V.Variable(j).N_etiquetas();
      if (nr1[j]<umbral){
	for (int s=0; s<longitud; s++)
	  cadena[tama+s]='1';
      }
      tama+=longitud;
    }


    double aux2;
    tama=0;
    for (int j=0; j<r1-1; j++)
    {
      aux2=nr1[j];
      longitud = V.Variable(j).N_etiquetas();

        if (nr1[j]>=umbral)
        {
	  for (int s=0; s<longitud; s++){
	    if (nb1[tama+s]=='0' and patr[tama+s]=='?') {
	      nb1[tama+s]='1';
	      cadena[tama+s]='1';
	      Bondad(G,V,E,R,valor2,N_valor,elemento_i,discriminante, reduccion_reglas, ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
		    basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar, !inGA, initial_rule, tam_rule_not_removed);
	      if (debug)
		cout << cadena << " v[0]=" << valor2[0] << " V[1]= " << valor2[1] << " V[2]= " << valor2[2] << " V[3]= " << valor2[3] <<  endl;
		  if ( (valor2[0]>valor[0] /*and valor2[0]>0*/) or
		       (valor2[0]==valor[0] and valor2[1]>valor[1]) or
		       (valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]>=valor[2]) /*or
		       (valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]==valor[2] and valor2[3]>=valor[3]) or
		       (valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]==valor[2] and valor2[3]==valor[3] and valor2[4]>=valor[4])*/
		    )
		  {
		  valor[0]=valor2[0];
		  valor[1]=valor2[1];
		  valor[2]=valor2[2];
		  valor[3]=valor2[3];
		  valor[4]=valor2[4];
		  //G.Valoracion(elemento_i,valor2);
		  mejora_global=true;
		  //cout << "  Anulado caso " << s << " de la variable " << V.Variable(j).Name() << "  Valor: " << valor2[0] << endl;
		  //cout << "\t";
		  if (debug){
		    cout << "\t\t\t\t" << cadena << endl;
		    cout << "  Anulado caso " << s << " de la variable " << V.Variable(j).Name() << "  Valor: " << valor2[0] << endl;
		    G.Valoracion(elemento_i,valor2);
		    G.PintaFitness_Estacionario(elemento_i);
		  }
	      }
	      else{
		  //cout << "Valor de caso " << s << " de la variable " << V.Variable(j).Name() << ": " << nb1[tama+s] << endl;
		  nb1[tama+s]='0';
		  cadena[tama+s]='0';
		  nr1[j]=aux2;
		  //Bondad(G,V,E,R,valor2,N_valor,elemento_i,discriminante, true, ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
		    //basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar);
	      //cout << "  Sigue a 0 el caso " << s << " de la variable " << V.Variable(j).Name() << "  Valor: " << valor2[0] << endl;
	      //cout << cadena << endl;
	      //cout << nb1 << endl;
	      }
	    }
	  }

	}
      tama+=longitud;
    }

    Bondad(G,V,E,R,valor,N_valor,elemento_i,discriminante, reduccion_reglas, ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
	   basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar,!inGA, initial_rule, tam_rule_not_removed);

    G.Valoracion(elemento_i,valor);

    if (debug){
      cout << "Valor Final\n";
      G.PintaFitness_Estacionario(elemento_i);
    }


    if (debug)
      cout << cadena << endl;


    /*cout << "Valor despues de Discriminacion: ";
    G.PintaFitness_Estacionario(elemento_i);


    cout << "cota de v[1]= " << nr1[1] << endl;
    nb1=G.Subpoblacion_Binaria(0,elemento_i,b1);
    nr1=G.Subpoblacion_Real(0,elemento_i,r1);
    cout << "cota de v[1]= " << nr1[1] << endl;

    cadena = nb1;
    tama=0;
    for (int j=0; j<r1-1; j++)
    {
      longitud = V.Variable(j).N_etiquetas();
      if (nr1[j]<umbral){
	for (int s=0; s<longitud; s++)
	  cadena[tama+s]='-';
      }
      tama+=longitud;
    }
    cout << "Resultado final: " << endl;
    //cout << cadena << endl;
    //cout << nb1 << endl;
    G.PintaFitness_Estacionario(elemento_i);*/

    /*if (mejora_global)
      cout << "\t\tMejorada la discriminacion de la regla\n";
    else
      cout << "\t\tNo hay mejora en la discriminacion de la regla\n";*/


    //char ch;
    //cin >> ch;

    delete [] valor;
    delete [] valor2;
    delete [] se_puede_anular;
}





void Obtener_Regla_Discriminante_E(multipoblacion &G, const VectorVar &V, example_set &E, ruleset &R, int elemento_i,
                               bool discriminante, bool reduccion_reglas, int* ponderacion_por_clase, double ***Tabla,
                               double ***Tabla_Nuevas, VectorVar &F, int *basica_para_regla,
			       int *errores_para_regla, int *reglas_a_eliminar, int &n_reglas_a_eliminar, bool debug, string patr,
			       bool inGA, string initial_rule, double nmas, double nmenos)
{
    int r1;
    double *nr1=G.Subpoblacion_Real(0,elemento_i,r1);
    double umbral = nr1[r1-1], aux;
    int N_valor=G.N_Valoracion();
    double *valor = new double [N_valor], *valor2 = new double[N_valor];
    string cadena;
    bool mejora=false, mejora_global=false;
    int tama=0;
    int b1;
    char *nb1=G.Subpoblacion_Binaria(0,elemento_i,b1);
    char *s= new char[b1+1];
    for (int i=0; i<b1; i++)
      s[i]=nb1[i];
    s[b1]='\0';
    cadena = s;
    delete [] s;
    int longitud;
    for (int j=0; j<r1-1; j++)
    {
      longitud = V.Variable(j).N_etiquetas();
      if (nr1[j]<umbral){
	for (int s=0; s<longitud; s++){
	  cadena[tama+s]='1';
	  nb1[tama+s]='1';
	}
      }
      tama+=longitud;
    }

    //Aciertos(R,V,E,F);

    Bondad_E(G,V,E,R,valor,N_valor,elemento_i,discriminante, reduccion_reglas, ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
          basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar,!inGA, initial_rule, nmas, nmenos);


    G.Valoracion(elemento_i,valor);
    //cout << cadena << endl;
    //cout << "Valor antes de Discriminacion:   ";
    if (debug){
      cout << "Valor antes de Discriminacion:   " << endl;
      cout << cadena << endl;
      cout << "Valor de partida\n";
      G.PintaFitness_Estacionario(elemento_i);
    }

    /*for (int i=0; i<r1; i++)
      cout << "\t" << nr1[i];
    cout << endl;*/
    //cout << "Umbral: " << umbral << endl;

    /*if (umbral<=0){
      umbral==0.01;
      nr1[r1-1]=umbral;
    }*/

    tama=0;

    bool *se_puede_anular = new bool[r1-1];

    for (int j=0; j<r1-1; j++)
    {
      longitud = V.Variable(j).N_etiquetas();
      if (nr1[j]>=umbral)
	se_puede_anular[j]=true;
      else
	se_puede_anular[j]=false;

      for (int s=0; s<longitud and se_puede_anular[j]; s++)
	if (patr[tama+s]=='0')
	  se_puede_anular[j]=false;

      if (patr[tama]=='x')
	nr1[j]=umbral-0.001;

      tama+=longitud;
    }

    for(int j=0; j<patr.length(); j++)
      if (patr[j]=='0' or patr[j]=='1')
	nb1[j]=patr[j];


    tama=0;
    for (int j=0; j<r1-1; j++)
    {
      longitud = V.Variable(j).N_etiquetas();

        if (nr1[j]>=umbral and se_puede_anular[j])
        {
            aux = nr1[j];
	    mejora=false;
            nr1[j]=umbral-0.001;

	    Bondad_E(G,V,E,R,valor2,N_valor,elemento_i,discriminante, reduccion_reglas, ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
		  basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar,!inGA, initial_rule, nmas, nmenos);
	    //cout << "\tValor[0]= " << valor2[0] << endl;
		  if ( (valor2[0]>valor[0] /*and valor2[0]>0*/) or
		       (valor2[0]==valor[0] and valor2[1]>valor[1]) or
		       (valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]>valor[2]) /*or
		       (valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]==valor[2] and valor2[3]>=valor[3]) or
		       (valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]==valor[2] and valor2[3]==valor[3] and valor2[4]>=valor[4])*/
		    )
		  {
		  valor[0]=valor2[0];
		  valor[1]=valor2[1];
		  valor[2]=valor2[2];
		  valor[3]=valor2[3];
		  valor[4]=valor2[4];
		  mejora=true;
		  mejora_global=true;
		  //G.Valoracion(elemento_i,valor2);
		  //cout << "\t\t\t" << nb1 << endl;
		  //cout << "Anulada la variable " << V.Variable(j).Name() << "  Valor: " << valor2[0] << endl;
		  if (debug){
		    //cout << "\t\t\t" << nb1 << endl;
		    cout << "Anulada la variable " << V.Variable(j).Name() << "  Valor: " << valor2[0] << endl;
		    G.Valoracion(elemento_i,valor2);
		    G.PintaFitness_Estacionario(elemento_i);
		  }
	    }
	    else{
	      nr1[j]=aux;
	      //Bondad(G,V,E,R,valor2,N_valor,elemento_i,discriminante, true, ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
		  //basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar);
	      //cout << "No se puede anular la variable " << V.Variable(j).Name() << "  Valor: " << valor2[0] << endl;
	    }

        }


      tama+=longitud;
    }

    tama=0;
    for (int j=0; j<r1-1; j++)
    {
      longitud = V.Variable(j).N_etiquetas();
      if (nr1[j]<umbral){
	for (int s=0; s<longitud; s++)
	  cadena[tama+s]='1';
      }
      tama+=longitud;
    }


    double aux2;
    tama=0;
    for (int j=0; j<r1-1; j++)
    {
      aux2=nr1[j];
      longitud = V.Variable(j).N_etiquetas();

        if (nr1[j]>=umbral)
        {
	  for (int s=0; s<longitud; s++){
	    if (nb1[tama+s]=='0' and patr[tama+s]=='?') {
	      nb1[tama+s]='1';
	      cadena[tama+s]='1';
	      Bondad_E(G,V,E,R,valor2,N_valor,elemento_i,discriminante, reduccion_reglas, ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
		    basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar, !inGA, initial_rule, nmas, nmenos);
	      if (debug)
		cout << cadena << " v[0]=" << valor2[0] << " V[1]= " << valor2[1] << " V[2]= " << valor2[2] << " V[3]= " << valor2[3] <<  endl;
		  if ( (valor2[0]>valor[0] /*and valor2[0]>0*/) or
		       (valor2[0]==valor[0] and valor2[1]>valor[1]) or
		       (valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]>=valor[2]) /*or
		       (valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]==valor[2] and valor2[3]>=valor[3]) or
		       (valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]==valor[2] and valor2[3]==valor[3] and valor2[4]>=valor[4])*/
		    )
		  {
		  valor[0]=valor2[0];
		  valor[1]=valor2[1];
		  valor[2]=valor2[2];
		  valor[3]=valor2[3];
		  valor[4]=valor2[4];
		  //G.Valoracion(elemento_i,valor2);
		  mejora_global=true;
		  //cout << "  Anulado caso " << s << " de la variable " << V.Variable(j).Name() << "  Valor: " << valor2[0] << endl;
		  //cout << "\t";
		  if (debug){
		    cout << "\t\t\t\t" << cadena << endl;
		    cout << "  Anulado caso " << s << " de la variable " << V.Variable(j).Name() << "  Valor: " << valor2[0] << endl;
		    G.Valoracion(elemento_i,valor2);
		    G.PintaFitness_Estacionario(elemento_i);
		  }
	      }
	      else{
		  //cout << "Valor de caso " << s << " de la variable " << V.Variable(j).Name() << ": " << nb1[tama+s] << endl;
		  nb1[tama+s]='0';
		  cadena[tama+s]='0';
		  nr1[j]=aux2;
		  //Bondad(G,V,E,R,valor2,N_valor,elemento_i,discriminante, true, ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
		    //basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar);
	      //cout << "  Sigue a 0 el caso " << s << " de la variable " << V.Variable(j).Name() << "  Valor: " << valor2[0] << endl;
	      //cout << cadena << endl;
	      //cout << nb1 << endl;
	      }
	    }
	  }

	}
      tama+=longitud;
    }

    Bondad_E(G,V,E,R,valor,N_valor,elemento_i,discriminante, reduccion_reglas, ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
	   basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar,!inGA, initial_rule, nmas, nmenos);

    G.Valoracion(elemento_i,valor);

    if (debug){
      cout << "Valor Final\n";
      G.PintaFitness_Estacionario(elemento_i);
    }


    if (debug)
      cout << cadena << endl;


    /*cout << "Valor despues de Discriminacion: ";
    G.PintaFitness_Estacionario(elemento_i);


    cout << "cota de v[1]= " << nr1[1] << endl;
    nb1=G.Subpoblacion_Binaria(0,elemento_i,b1);
    nr1=G.Subpoblacion_Real(0,elemento_i,r1);
    cout << "cota de v[1]= " << nr1[1] << endl;

    cadena = nb1;
    tama=0;
    for (int j=0; j<r1-1; j++)
    {
      longitud = V.Variable(j).N_etiquetas();
      if (nr1[j]<umbral){
	for (int s=0; s<longitud; s++)
	  cadena[tama+s]='-';
      }
      tama+=longitud;
    }
    cout << "Resultado final: " << endl;
    //cout << cadena << endl;
    //cout << nb1 << endl;
    G.PintaFitness_Estacionario(elemento_i);*/

    /*if (mejora_global)
      cout << "\t\tMejorada la discriminacion de la regla\n";
    else
      cout << "\t\tNo hay mejora en la discriminacion de la regla\n";*/


    //char ch;
    //cin >> ch;

    delete [] valor;
    delete [] valor2;
    delete [] se_puede_anular;
}






string Patron1Regla(ruleset &R_1, const VectorVar &V, int pos){
  genetcode r_1, r_2, *aux1, *aux2;
  aux1=R_1.Get_Code_Rule(pos);
  r_1=*aux1;

 int b,*nb,n,*nn,**nnn,r,*nr;
 char **mb;
 double **mr;
 string regla1;

 r_1.GetBinary(b,nb,mb);
 r_1.GetReal(r,nr,mr);

 char *s= new char[nb[0]+1];
 for (int j=0; j<nb[0]; j++)
   s[j]=mb[0][j];
 s[nb[0]]='\0';
 regla1=s;
 delete []s;

 int tama=0, longitud;

 /*for (int j=0;j<V.N_Antecedente(); j++){
     longitud = V.SizeDomain(j);
     if (mr[0][j]<mr[0][nr[0]-1]){
       for (int s=0; s<longitud; s++)
	 regla1[tama+s]='1';
     }
   tama+=longitud;
 }*/

 return regla1;

}


string Patron2Patrones(string pat1, string pat2){
  string pat3=pat1;

 for (int i=0; i<pat1.length(); i++){
   if (pat1[i]==pat2[i])
     pat3[i]=pat1[i];
   else if  ( (pat1[i]=='1' and pat2[i]=='?') or (pat1[i]=='?' and pat2[i]=='1') )
     pat3[i]='1';
   else if  ( (pat1[i]=='0' and pat2[i]=='?') or (pat1[i]=='?' and pat2[i]=='0') )
     pat3[i]='0';
   else
     pat3[i]='?';
 }

 return pat3;

}


string Patron2Reglas(ruleset &R_1, ruleset &R_2, const VectorVar &V, int pos){
  genetcode r_1, r_2, *aux1, *aux2;
  aux1=R_1.Get_Code_Rule(pos);
  r_1=*aux1;
  aux2=R_2.Get_Code_Rule(pos);
  r_2=*aux2;

 int b,*nb,n,*nn,**nnn,r,*nr;
 char **mb;
 double **mr;
 string regla1, regla2, regla3;

 r_1.GetBinary(b,nb,mb);
 r_1.GetReal(r,nr,mr);

 char *s= new char[nb[0]+1];
 for (int j=0; j<nb[0]; j++)
   s[j]=mb[0][j];
 s[nb[0]]='\0';
 regla1=s;
 delete []s;

 //cout << regla1 << endl;

 int longitud, tama=0;

 for (int j=0;j<V.N_Antecedente(); j++){
     longitud = V.SizeDomain(j);
     if (mr[0][j]<mr[0][nr[0]-1]){
       for (int s=0; s<longitud; s++)
	 regla1[tama+s]='1';
     }
   tama+=longitud;
 }

 //cout << regla1 << endl;

 r_2.GetBinary(b,nb,mb);
 r_2.GetReal(r,nr,mr);

 s= new char[nb[0]+1];
 for (int j=0; j<nb[0]; j++)
   s[j]=mb[0][j];
 s[nb[0]]='\0';
 regla2=s;
 delete []s;

 //cout << regla2 << endl;

 tama=0;

 for (int j=0;j<V.N_Antecedente(); j++){
     longitud = V.SizeDomain(j);
     if (mr[0][j]<mr[0][nr[0]-1]){
       for (int s=0; s<longitud; s++)
	 regla2[tama+s]='1';
     }
   tama+=longitud;
 }

 //cout << regla2 << endl;

 regla3=regla1;
 for (int i=0; i<nb[0]; i++){
   if (regla1[i]==regla2[i])
     regla3[i]=regla1[i];
   else
     regla3[i]='?';
 }

 return regla3;

}

void Pinta2Reglas(ruleset &R_1, ruleset &R_2, const VectorVar &V, int pos){
  genetcode r_1, r_2, *aux;
  aux=R_1.Get_Code_Rule(pos);
  r_1=*aux;
  aux=R_2.Get_Code_Rule(pos);
  r_2=*aux;

 int b,*nb,n,*nn,**nnn,r,*nr;
 char **mb;
 double **mr;
 string regla1, regla2, regla3;

 r_1.GetBinary(b,nb,mb);
 r_1.GetReal(r,nr,mr);

 char *s= new char[nb[0]+1];
 for (int j=0; j<nb[0]; j++)
   s[j]=mb[0][j];
 s[nb[0]]='\0';
 regla1=s;
 delete []s;

 int longitud, tama=0;

 for (int j=0;j<V.N_Antecedente(); j++){
     longitud = V.SizeDomain(j);
     if (mr[0][j]<mr[0][nr[0]-1]){
       for (int s=0; s<longitud; s++)
	 regla1[tama+s]='-';
     }
   tama+=longitud;
 }


 r_2.GetBinary(b,nb,mb);
 r_2.GetReal(r,nr,mr);

 s= new char[nb[0]+1];
 for (int j=0; j<nb[0]; j++)
   s[j]=mb[0][j];
 s[nb[0]]='\0';
 regla2=s;
 delete []s;

 tama=0;

 for (int j=0;j<V.N_Antecedente(); j++){
     longitud = V.SizeDomain(j);
     if (mr[0][j]<mr[0][nr[0]-1]){
       for (int s=0; s<longitud; s++)
	 regla2[tama+s]='-';
     }
   tama+=longitud;
 }


 cout << regla1 << endl;
 cout << regla2 << endl;

 regla3=regla1;
 for (int i=0; i<nb[0]; i++){
   if (regla1[i]== regla2[i] and regla1[i]=='-')
     regla3[i]='x';
   else if (regla1[i]==regla2[i])
     regla3[i]=regla1[i];
   else if (regla1[i]=='1' and regla2[i]=='-')
     regla3[i]=regla1[i];
   else
     regla3[i]='?';

 }

 cout << regla3 << endl;


 //char ch;
 //cin >> ch;

}




void Obtener_Regla_Descriptiva(multipoblacion &G, const VectorVar &V, example_set &E, ruleset &R, int elemento_i,
                               bool discriminante, bool reduccion_reglas, int* ponderacion_por_clase, double ***Tabla,
                               double ***Tabla_Nuevas, VectorVar &F, int *basica_para_regla,
			       int *errores_para_regla, int *reglas_a_eliminar, int &n_reglas_a_eliminar, bool debug, string patr,
			       bool inGA, string initial_rule, int tam_rule_not_removed)
{
    int r1;
    double *nr1=G.Subpoblacion_Real(0,elemento_i,r1);
    double umbral = nr1[r1-1], aux;
    int N_valor=G.N_Valoracion();
    double *valor = new double [N_valor], *valor2 = new double[N_valor];
    bool mejora, mejora_global=false;
    int tama=0;
    int b1;
    char *nb1=G.Subpoblacion_Binaria(0,elemento_i,b1);


    char *s= new char[b1+1];
    for (int i=0; i<b1; i++)
      s[i]=nb1[i];
    s[b1]='\0';
    string cadena=s;
    delete [] s;
    int longitud;

    for (int j=0; j<r1-1; j++)
    {
      longitud = V.Variable(j).N_etiquetas();
      if (nr1[j]<umbral){
	for (int s=0; s<longitud; s++)
	  cadena[tama+s]='1';
      }
      tama+=longitud;
    }

    /*for (int i=0; i<r1; i++)
      cout << "\t" << nr1[i];
    cout << endl;*/
    //cout << "Umbral: " << umbral << endl;

    //Aciertos(R,V,E,F);

    Bondad(G,V,E,R,valor,N_valor,elemento_i,discriminante, reduccion_reglas, ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
          basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar, !inGA, initial_rule, tam_rule_not_removed);

    //if (!inGA and ( !valor[0]==0 or !valor[1]==0) ){

      G.Valoracion(elemento_i,valor);
      if (debug){
	cout << "Valor al principio\n";
	G.PintaFitness_Estacionario(elemento_i);
      }

      //G.PintaFitness_Estacionario(elemento_i);
      //cout << cadena << endl;

      //for(int j=0; j<patr.length(); j++)
      //  nb1[j]='1';

      //cout << "Patron: " << patr << endl;

      for(int j=0; j<patr.length(); j++)
	if (patr[j]=='0' or patr[j]=='1')
	  nb1[j]=patr[j];


      tama=0;
      for (int j=0; j<r1-1; j++)
      {
	longitud = V.Variable(j).N_etiquetas();

	  if (nr1[j]<umbral and patr[tama]!='x')
	  {
	      aux = nr1[j];
	      mejora=false;
	      nr1[j]=umbral;
	      for (int s=0; s<longitud; s++){
		if (patr[tama+s]=='?'){
		  nb1[tama+s]='1';
		  cadena[tama+s]='1';
		}
		else {
		  nb1[tama+s]=patr[tama+s];
		  cadena[tama+s]=patr[tama+s];
		}
	      }

	      for (int s=0; s<longitud; s++){
		if (patr[tama+s]=='?' and nb1[tama+s]=='1'){
		    nb1[tama+s]='0';
		    cadena[tama+s]='0';
		    Bondad(G,V,E,R,valor2,N_valor,elemento_i,discriminante, reduccion_reglas, ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
			  basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar,!inGA, initial_rule, tam_rule_not_removed);
		    if (debug)
		      cout << cadena << " v[0]=" << valor2[0] << " V[1]= " << valor2[1] << " V[2]= " << valor2[2] << " V[3]= " << valor2[3] <<  endl;

		    if ( (valor2[0]>valor[0] /*and valor2[0]>0*/) or
			(valor2[0]==valor[0] and valor2[1]>valor[1]) or
			(valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]>=valor[2]) /*or
			(valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]==valor[2] and valor2[3]>=valor[3]) or
			(valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]==valor[2] and valor2[3]==valor[3] and valor2[4]>=valor[4])*/
		      )
		    {
		      valor[0]=valor2[0];
		      valor[1]=valor2[1];
		      valor[2]=valor2[2];
		      valor[3]=valor2[3];
		      valor[4]=valor2[4];
		      mejora=true;
		      mejora_global=true;
		      G.Valoracion(elemento_i,valor2);
		      if (debug){
			G.PintaFitness_Estacionario(elemento_i);
		      }

		    //G.PintaFitness_Estacionario(elemento_i);
		    //cout << cadena << endl;
		    }
		    else{
			nb1[tama+s]='1';
			cadena[tama+s]='1';
			nr1[j]=umbral;
		    }
		}
	      }

	      if (!mejora){
		nr1[j]=aux;
		for (int s=0; s<longitud; s++){
		  cadena[tama+s]='-';
		}
	      }
	  }
	  else if (patr[tama]!='x') {
	      for (int s=0; s<longitud; s++)
		if (nb1[tama+s]=='1' and patr[tama+s]=='?') {
		  nb1[tama+s]='0';
		  cadena[tama+s]='0';
		  Bondad(G,V,E,R,valor2,N_valor,elemento_i,discriminante, reduccion_reglas, ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
			basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar,!inGA, initial_rule, tam_rule_not_removed);
		    if ( (valor2[0]>valor[0] /*and valor2[0]>0*/) or
			(valor2[0]==valor[0] and valor2[1]>valor[1]) or
			(valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]>valor[2]) /* or
			(valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]==valor[2] and valor2[3]>=valor[3]) or
			(valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]==valor[2] and valor2[3]==valor[3] and valor2[4]>=valor[4])*/
		      )
		    {
		      valor[0]=valor2[0];
		      valor[1]=valor2[1];
		      valor[2]=valor2[2];
		      valor[3]=valor2[3];
		      valor[4]=valor2[4];
		      mejora=true;
		      mejora_global=true;
		      G.Valoracion(elemento_i,valor2);
		      if (debug){
			G.PintaFitness_Estacionario(elemento_i);
		    }

		      //G.PintaFitness_Estacionario(elemento_i);
		      //cout << cadena << endl;
		  }
		  else{
		      nb1[tama+s]='1';
		      cadena[tama+s]='1';
		  }
		}
	  }

	tama+=longitud;
      }

      Bondad(G,V,E,R,valor,N_valor,elemento_i,discriminante, reduccion_reglas, ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
	    basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar, !inGA, initial_rule, tam_rule_not_removed);

      G.Valoracion(elemento_i,valor);

      if (debug){
	cout << "Valor final\n";
	G.PintaFitness_Estacionario(elemento_i);
      }

      /*cout << "Resultado final: " << endl;
      nb1=G.Subpoblacion_Binaria(0,elemento_i,b1);
      cadena=nb1;
      tama=0;
      for (int j=0; j<r1-1; j++)
      {
	longitud = V.Variable(j).N_etiquetas();
	if (nr1[j]<umbral){
	  for (int s=0; s<longitud; s++)
	    cadena[tama+s]='-';
	}
	tama+=longitud;
      }
      cout << cadena << endl;
      G.PintaFitness_Estacionario(elemento_i);    */

      //char ch;
      //cin >> ch;
      /*if (mejora_global)
	cout << "\t\tMejorada la descripcion de la regla\n";
      else
	cout << "\t\tNo hay mejora en la descripcion de la regla\n";*/
    //}

    delete [] valor;
    delete [] valor2;
}




void Obtener_Regla_Descriptiva_E(multipoblacion &G, const VectorVar &V, example_set &E, ruleset &R, int elemento_i,
                               bool discriminante, bool reduccion_reglas, int* ponderacion_por_clase, double ***Tabla,
                               double ***Tabla_Nuevas, VectorVar &F, int *basica_para_regla,
			       int *errores_para_regla, int *reglas_a_eliminar, int &n_reglas_a_eliminar, bool debug, string patr,
			       bool inGA, string initial_rule, double nmas, double nmenos)
{
    int r1;
    double *nr1=G.Subpoblacion_Real(0,elemento_i,r1);
    double umbral = nr1[r1-1], aux;
    int N_valor=G.N_Valoracion();
    double *valor = new double [N_valor], *valor2 = new double[N_valor];
    bool mejora, mejora_global=false;
    int tama=0;
    int b1;
    char *nb1=G.Subpoblacion_Binaria(0,elemento_i,b1);


    char *s= new char[b1+1];
    for (int i=0; i<b1; i++)
      s[i]=nb1[i];
    s[b1]='\0';
    string cadena=s;
    delete [] s;
    int longitud;

    for (int j=0; j<r1-1; j++)
    {
      longitud = V.Variable(j).N_etiquetas();
      if (nr1[j]<umbral){
	for (int s=0; s<longitud; s++)
	  cadena[tama+s]='1';
      }
      tama+=longitud;
    }

    /*for (int i=0; i<r1; i++)
      cout << "\t" << nr1[i];
    cout << endl;*/
    //cout << "Umbral: " << umbral << endl;

    //Aciertos(R,V,E,F);

    Bondad_E(G,V,E,R,valor,N_valor,elemento_i,discriminante, reduccion_reglas, ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
          basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar, !inGA, initial_rule, nmas, nmenos);

    //if (!inGA and ( !valor[0]==0 or !valor[1]==0) ){

      G.Valoracion(elemento_i,valor);
      if (debug){
	cout << "Valor al principio\n";
	G.PintaFitness_Estacionario(elemento_i);
      }

      //G.PintaFitness_Estacionario(elemento_i);
      //cout << cadena << endl;

      //for(int j=0; j<patr.length(); j++)
      //  nb1[j]='1';

      //cout << "Patron: " << patr << endl;

      for(int j=0; j<patr.length(); j++)
	if (patr[j]=='0' or patr[j]=='1')
	  nb1[j]=patr[j];


      tama=0;
      for (int j=0; j<r1-1; j++)
      {
	longitud = V.Variable(j).N_etiquetas();

	  if (nr1[j]<umbral and patr[tama]!='x')
	  {
	      aux = nr1[j];
	      mejora=false;
	      nr1[j]=umbral;
	      for (int s=0; s<longitud; s++){
		if (patr[tama+s]=='?'){
		  nb1[tama+s]='1';
		  cadena[tama+s]='1';
		}
		else {
		  nb1[tama+s]=patr[tama+s];
		  cadena[tama+s]=patr[tama+s];
		}
	      }

	      for (int s=0; s<longitud; s++){
		if (patr[tama+s]=='?' and nb1[tama+s]=='1'){
		    nb1[tama+s]='0';
		    cadena[tama+s]='0';
		    Bondad_E(G,V,E,R,valor2,N_valor,elemento_i,discriminante, reduccion_reglas, ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
			  basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar,!inGA, initial_rule, nmas, nmenos);
		    if (debug)
		      cout << cadena << " v[0]=" << valor2[0] << " V[1]= " << valor2[1] << " V[2]= " << valor2[2] << " V[3]= " << valor2[3] <<  endl;

		    if ( (valor2[0]>valor[0] /*and valor2[0]>0*/) or
			(valor2[0]==valor[0] and valor2[1]>valor[1]) or
			(valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]>=valor[2]) /*or
			(valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]==valor[2] and valor2[3]>=valor[3]) or
			(valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]==valor[2] and valor2[3]==valor[3] and valor2[4]>=valor[4])*/
		      )
		    {
		      valor[0]=valor2[0];
		      valor[1]=valor2[1];
		      valor[2]=valor2[2];
		      valor[3]=valor2[3];
		      valor[4]=valor2[4];
		      mejora=true;
		      mejora_global=true;
		      G.Valoracion(elemento_i,valor2);
		      if (debug){
			G.PintaFitness_Estacionario(elemento_i);
		      }

		    //G.PintaFitness_Estacionario(elemento_i);
		    //cout << cadena << endl;
		    }
		    else{
			nb1[tama+s]='1';
			cadena[tama+s]='1';
			nr1[j]=umbral;
		    }
		}
	      }

	      if (!mejora){
		nr1[j]=aux;
		for (int s=0; s<longitud; s++){
		  cadena[tama+s]='-';
		}
	      }
	  }
	  else if (patr[tama]!='x') {
	      for (int s=0; s<longitud; s++)
		if (nb1[tama+s]=='1' and patr[tama+s]=='?') {
		  nb1[tama+s]='0';
		  cadena[tama+s]='0';
		  Bondad_E(G,V,E,R,valor2,N_valor,elemento_i,discriminante, reduccion_reglas, ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
			basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar,!inGA, initial_rule, nmas, nmenos);
		    if ( (valor2[0]>valor[0] /*and valor2[0]>0*/) or
			(valor2[0]==valor[0] and valor2[1]>valor[1]) or
			(valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]>valor[2]) /* or
			(valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]==valor[2] and valor2[3]>=valor[3]) or
			(valor2[0]==valor[0] and valor2[1]==valor[1] and valor2[2]==valor[2] and valor2[3]==valor[3] and valor2[4]>=valor[4])*/
		      )
		    {
		      valor[0]=valor2[0];
		      valor[1]=valor2[1];
		      valor[2]=valor2[2];
		      valor[3]=valor2[3];
		      valor[4]=valor2[4];
		      mejora=true;
		      mejora_global=true;
		      G.Valoracion(elemento_i,valor2);
		      if (debug){
			G.PintaFitness_Estacionario(elemento_i);
		    }

		      //G.PintaFitness_Estacionario(elemento_i);
		      //cout << cadena << endl;
		  }
		  else{
		      nb1[tama+s]='1';
		      cadena[tama+s]='1';
		  }
		}
	  }

	tama+=longitud;
      }

      Bondad_E(G,V,E,R,valor,N_valor,elemento_i,discriminante, reduccion_reglas, ponderacion_por_clase, R.N_rule(), Tabla, Tabla_Nuevas, F,
	    basica_para_regla, errores_para_regla, reglas_a_eliminar, n_reglas_a_eliminar, !inGA, initial_rule, nmas, nmenos);

      G.Valoracion(elemento_i,valor);

      if (debug){
	cout << "Valor final\n";
	G.PintaFitness_Estacionario(elemento_i);
      }

      /*cout << "Resultado final: " << endl;
      nb1=G.Subpoblacion_Binaria(0,elemento_i,b1);
      cadena=nb1;
      tama=0;
      for (int j=0; j<r1-1; j++)
      {
	longitud = V.Variable(j).N_etiquetas();
	if (nr1[j]<umbral){
	  for (int s=0; s<longitud; s++)
	    cadena[tama+s]='-';
	}
	tama+=longitud;
      }
      cout << cadena << endl;
      G.PintaFitness_Estacionario(elemento_i);    */

      //char ch;
      //cin >> ch;
      /*if (mejora_global)
	cout << "\t\tMejorada la descripcion de la regla\n";
      else
	cout << "\t\tNo hay mejora en la descripcion de la regla\n";*/
    //}

    delete [] valor;
    delete [] valor2;
}



void Recalcular_Basicos (ruleset &R, example_set &E, double *adapt, double peso_r){
double **Rule = new double*[R.N_rule()];
double *peso = new double[R.N_rule()];
int aciertos=0, fallos=0, disparada, mejor_valor;

  for (int i=0; i<R.N_rule();i++){
    Rule[i]=R.SeeAdaptaciones_of_rule(i);
    peso[i]=R.Get_Weight(i);
  }

  for (int j=0; j<E.N_Examples(); j++){
    disparada=-1;
    mejor_valor=0;
    for (int i=0; i<R.N_rule();i++){
      if (Rule[i][j]>fabs(mejor_valor) or (Rule[i][j]==fabs(mejor_valor) and fabs(mejor_valor)>0 and peso[i]>peso[disparada]) ){
        disparada=i;
        mejor_valor=Rule[i][j];
      }
    }
  }

delete [] Rule;
delete [] peso;
}




double Diferencia_Adaptacion_Negativa(double *esta_regla, double *regla_aprendida, int tama){
  double contador=0;
  for (int i=0; i<tama; i++)
    if (esta_regla[i]<=0 and esta_regla[i]>regla_aprendida[i])
       contador+=(esta_regla[i]-regla_aprendida[i]);
    else if (regla_aprendida[i]<=0 and regla_aprendida[i]>esta_regla[i])
            contador+=(esta_regla[i]-regla_aprendida[i]);

  return contador;
}



double Weight_Individual(multipoblacion &G, const VectorVar &V, const example_set &E, VectorVar &F, double ***Tabla, double ***Tabla_Nuevas, int elemento_i,
			 double &n_mas, double &n_menos){

  int b1, b2, n1, n2, r1, r2, r3;
    char *nb1, *nb2;
    int *nn1,*nn2,r;
    double *nr1, *nr2, *nr3;
    string regla, regla_funciones;
    int var1, var2;
    int ne=E.N_Examples();
    double peso;
    int j;

    nb1=G.Subpoblacion_Binaria(0,elemento_i,b1);
    char *s= new char[b1+1];
    for (j=0; j<b1; j++)
    {
        s[j]=nb1[j];
        s[b1]='\0';
    }
    regla=s;
    delete []s;


    nb2=G.Subpoblacion_Binaria(1,elemento_i,b2);
    s= new char[b2+1];
    for (j=0; j<b2; j++)
    {
        s[j]=nb2[j];
        s[b2]='\0';
    }
    regla_funciones=s;
    delete []s;


    nn1=G.Subpoblacion_Entera(0,elemento_i,n1);
    nn2=G.Subpoblacion_Entera(1,elemento_i,n2);
    nr1=G.Subpoblacion_Real(0,elemento_i,r1);

    nr2 = new double[ne];

  double positivos=0, negativos=0, aux_p, aux_n;

  for (j=0; j<E.N_Examples(); j++)
    {

	T_AdaptacionC(V,Tabla,j,nn1[0],aux_p,aux_n);

	//esta_cubierto = E.Is_Covered(j);

	double x = T_Adaptacion_Antecedente(V,Tabla,j,regla,nr1,nr1[r1-1]);

	nr2[j]=x;
	//nr2[j]=1;
	nr2[j]=T_Adaptacion_Antecedente_Nuevas(F, Tabla_Nuevas, j, regla_funciones, nn2, n2, nr2[j]);


	positivos=positivos+(nr2[j]*aux_p);
	negativos=negativos+(nr2[j]*aux_n);
	if (aux_p<aux_n)
	    nr2[j]=-nr2[j];

    }

  n_mas = positivos;
  n_menos = negativos;

  peso = (positivos+1) / (positivos+negativos+1);

  delete [] nr2;

  return peso;
}



double Bondad (multipoblacion &G, const VectorVar &V, const example_set &E, ruleset &R, double *valor, int N_valor,
               int elemento_i, bool discriminante, bool reducir_reglas, int *ponderacion_por_clase, int nrule,
               double ***Tabla, double ***Tabla_Nuevas, VectorVar &F, int *basica_para_regla, int *errores_para_regla, int *reglas_a_eliminar, int &n_reglas_a_eliminar,
	       bool esp_gen, string initial_rule, int tam_rule_not_removed)
{

    int b1, b2, n1, n2, r1, r2, r3;
    char *nb1, *nb2;
    int *nn1,*nn2,r;
    double *nr1, *nr2, *nr3;
    string regla, regla_funciones;
    int var1, var2;
    int ne=E.N_Examples();
    vectordouble w(V.N_Antecedente());
    double positivos=0, negativos=0, rango, eficacia, negativos_que_conllevan_error=0, positivos_a_incrementar=0,
                                  positivos_sobrecubrimiento=0, negativos2=0;
    int j;
    double aciertos=0, fallos=0, aciertos_aqui=0, fallos_aqui=0, sin_cubrir=0, aciertos_no_cubiertos=0, max_pos=0,
                              max_neg=0, max_neg2=0;
    int clase = G.ClaseIndividuo(elemento_i);
    bool menos_infinito=true, valida=true;
    int variables_irrelevantes=0, funciones_irrelevantes=0, variables_estables=0, funciones_estables=0,
        ceros_variables=0, ceros_funciones=0, unos_variables=0, unos_funciones=0,
	variables_relevantes, funciones_relevantes;

    nrule= R.N_rule();


    nb1=G.Subpoblacion_Binaria(0,elemento_i,b1);
    char *s= new char[b1+1];
    valor[2]=0;  //numero de unos en el cromosoma
    for (j=0; j<b1; j++)
    {
        s[j]=nb1[j];
        if (nb1[j]=='1')
            valor[2]++;
        s[b1]='\0';
    }
    regla=s;
    delete []s;


    nb2=G.Subpoblacion_Binaria(1,elemento_i,b2);
    s= new char[b2+1];
    for (j=0; j<b2; j++)
    {
        s[j]=nb2[j];
        if (nb2[j]=='1')
            valor[2]++;
        s[b2]='\0';
    }
    regla_funciones=s;
    delete []s;


    nn1=G.Subpoblacion_Entera(0,elemento_i,n1);
    nn2=G.Subpoblacion_Entera(1,elemento_i,n2);
    nr1=G.Subpoblacion_Real(0,elemento_i,r1);
    nr2 = new double[ne];


    valor[1]=0; //numero de variable irrelevantes
    valor[2]=0; //funciones irrelevantes
    valor[3]=0; //asignaciones comprensibles
    valor[4]=0; //numero de ceros
    int pos=0; // posicion en el cromosoma nivel de valor
    int tamanio; // numero de genes de una variable

    //cout << "Umbral: " << nr1[r1-1] << endl;
    for (j=0; j<r1-1; j++)
    {
//cout << j << endl;
        tamanio= V.Variable(j).N_etiquetas();
        if (nr1[j]<nr1[r1-1])
            variables_irrelevantes++;
        else
        {
           if (Relevante(&(nb1[pos]),tamanio)){
              if (Asignacion_Comprensible(&(nb1[pos]),tamanio)){
                variables_estables++;
                //ceros_variables+=NumeroCeros(&(nb1[pos]),tamanio);
                //unos_variables+=NumeroUnos(&(nb1[pos]),tamanio);
              }
              //else {
                ceros_variables+=NumeroCeros(&(nb1[pos]),tamanio);
                unos_variables+=NumeroUnos(&(nb1[pos]),tamanio);
                // Asignacion no comprensible
                // considero entonces la regla como invalida
                //valida=false;
              //}
           }
           else {
              nr1[j]=nr1[r1-1]-0.001;
              variables_irrelevantes++;
           }

        }
        pos=pos+tamanio;
    }

    // Quitar funciones repetidas

    for (int i=1; i<n2; i++){
      if (nn2[i]!=0) {
        int t=0;
        while (t<i and nn2[t]!=nn2[i])
          t++;
        if (t<i)
          nn2[i]=0;
      }
    }

    variables_relevantes = (r1-1) - variables_irrelevantes;


    pos=0;
    for (j=0; j<n2; j++)
    {
        if (nn2[j]!=0)
        {
           tamanio= F.SizeDomain(nn2[j]-1);
           if (Relevante(&(nb2[pos]),tamanio)){
              if (Asignacion_Comprensible(&(nb2[pos]),tamanio)){
                 funciones_estables++;
                 ceros_funciones+=NumeroCeros(&(nb2[pos]),tamanio);
                 unos_funciones+=NumeroUnos(&(nb2[pos]),tamanio);
              }
              else {
                // Asignacion no comprensible
                // considero entonces la regla como invalida
                //valida=false;
              }
           }
           else {
              funciones_irrelevantes++;
              nn2[j]=0;
           }
        }
        else {
            tamanio=2;
            funciones_irrelevantes++;
        }
        tamanio=2;
        pos=pos+tamanio;
    }

    funciones_relevantes = (n2) - funciones_irrelevantes;


//Distingue entre descripciones discriminantes o descriptivas
    /*if (!discriminante)
        valor[1]=r1-1+n2-valor[1]-valor[2];*/


    bool activada_memoria=false;
    int *basica_esta_regla = 0;
    int *errores_esta_regla = 0;
    if (nrule>0) {
      basica_esta_regla = new int[nrule];
      errores_esta_regla = new int[nrule];
      activada_memoria=true;
      for (int i=0; i<nrule; i++){
        basica_esta_regla[i]=0;
        errores_esta_regla[i]=0;
      }
    }

    int *ejemplos_afectado_por_esta_regla = new int[ne];

//double *aux = new double[ne];
    double aux_p;
    double aux_n;
    double peso, kk;
    double fitness0_global_regla=-ne;
// Comprobar que es una regla valida
    bool esta_cubierto;
//bool regla_valida = V.Es_Valida(regla,nr1,nr1[r1-1],valor[2]);
    bool regla_valida = V.Es_Valida(regla,nr1,nr1[r1-1],kk);
    double afecta=0;

    double n_mas=0, n_menos=0;

    if (valida and !esp_gen)
       valida = ponderacion_por_clase[clase]!=0;
    else {
       valida = true;
       if (!regla_valida){
	 //cout << "Regla considerada como no valida\n";
       }
    }

// Calculo previo para obtener el peso de la regla sobre atributos originales
    if (regla_valida and valida)
    {
        for (j=0; j<ne; j++)
        {

            T_AdaptacionC(V,Tabla,j,nn1[0],aux_p,aux_n);

            //esta_cubierto = E.Is_Covered(j);

            nr2[j]=T_Adaptacion_Antecedente(V,Tabla,j,regla,nr1,nr1[r1-1]);
            //nr2[j]=1;
            nr2[j]=T_Adaptacion_Antecedente_Nuevas(F, Tabla_Nuevas, j, regla_funciones, nn2, n2, nr2[j]);

            if (nr2[j]*aux_p>0)
              afecta++;

            positivos=positivos+(nr2[j]*aux_p);
            negativos=negativos+(nr2[j]*aux_n);
            if (aux_p<aux_n)
                nr2[j]=-nr2[j];

        }

        //for (int i=0; i<nrule; i++)
          //cout <<"Regla " << i <<": " << basica_para_regla[i] << endl;
        //fitness0_global_regla = (positivos_global-negativos_global)/(ponderacion_por_clase[clase]);
    }


//peso = (positivos+1) / (positivos+negativos+1);
//peso = 1.0;
    n_mas=positivos;
    n_menos=negativos;
    peso = (positivos+1) / (positivos+negativos+1);
//peso = nr3[0];

 if (positivos>0)
   menos_infinito=false;


 /*if (afecta< 0.20*ponderacion_por_clase[clase])
   valida=false;*/

 if (ponderacion_por_clase[clase]>0)
   afecta = afecta / ponderacion_por_clase[clase];

 double positivos_global=0, negativos_global=0;
// Dado el peso de la regla obtener la bondad de la misma.
    int positivos_aux=0;
    if (regla_valida and valida)
    {
        positivos=0;
        negativos=0;
        for (j=0; j<ne; j++)
        {
            if (nr2[j]>=0)
              positivos_global+= (nr2[j]*peso);
            else
              negativos_global+= (-nr2[j]*peso);

            ejemplos_afectado_por_esta_regla[j]=0;

            esta_cubierto = E.Is_Covered(j);
            if (!esta_cubierto)
            {
                if (nr2[j]>0)
                {
                    if ( (nr2[j]*peso>E.Grade_Is_Negative_Covered(j) ) or (nr2[j]*peso==E.Grade_Is_Negative_Covered(j) and peso>E.Peso_Is_Negative_Covered(j) ) )
                    {
                        positivos=positivos+(nr2[j]*peso);
                        aciertos++;
                        ejemplos_afectado_por_esta_regla[j]=1;
                        if (nrule>0 and E.Rule_Aislada(j)==-1 and E.Rule_Used(j)>=0 and E.Rule_Used(j)<nrule)
                          errores_esta_regla[E.Rule_Used(j)]++;
                    }
                    else {
                        positivos_a_incrementar += (nr2[j]*peso);
                    }
                }
                else if (nr2[j]<0) {

                       if (  (-nr2[j]*peso>E.Grade_Is_Positive_Covered(j)) or (-nr2[j]*peso==E.Grade_Is_Positive_Covered(j) and peso>E.Peso_Is_Positive_Covered(j)) ){
                          //negativos=negativos+(-nr2[j]*peso);
                          //if (E.Grade_Is_Positive_Covered(j)>0)
                            negativos=negativos+((-nr2[j]*peso)/ponderacion_por_clase[E.Class_of_example(V.Consecuente(),j)]);
                          //negativos=negativos+((-nr2[j]*peso)/ponderacion_por_clase[E.Class_of_example(V.Consecuente(),j)]);
                          ejemplos_afectado_por_esta_regla[j]=-1;

                          if (nrule>0 and E.Rule_Aislada(j)==-1 and E.Rule_Used(j)>=0 and E.Rule_Used(j)<nrule)
                            errores_esta_regla[E.Rule_Used(j)]++;
                       }

                }
            }
            else
            {
                if (nr2[j]>0)
                {
                   if (nr2[j]*peso>E.Grade_Is_Negative_Covered(j) or
                      (nr2[j]*peso==E.Grade_Is_Negative_Covered(j) and peso>E.Peso_Is_Negative_Covered(j) ) ){
                        ejemplos_afectado_por_esta_regla[j]=1;
                        if (nrule>0 and E.Rule_Aislada(j)==1 and E.Rule_Used(j)>=0 and E.Rule_Used(j)<nrule)
                           basica_esta_regla[E.Rule_Used(j)]++;
                        positivos_sobrecubrimiento+=1;
                   }

                }
                else  if (nr2[j]<0 and ( (-nr2[j])*peso>E.Grade_Is_Positive_Covered(j) or
                         ((-nr2[j])*peso==E.Grade_Is_Positive_Covered(j) and peso > E.Peso_Is_Positive_Covered(j) ) ) )
                {
                    ejemplos_afectado_por_esta_regla[j]=-1;
                    //negativos=negativos+(-nr2[j]*peso);
                    negativos=negativos+((-nr2[j]*peso)/ponderacion_por_clase[E.Class_of_example(V.Consecuente(),j)]);
                    fallos++;
                    if (nrule>0 and E.Rule_Aislada(j)==1 and E.Rule_Used(j)>=0 and E.Rule_Used(j)<nrule)
                      basica_esta_regla[E.Rule_Used(j)]++;
                }
            }
        }
    }
// La regla es invalida
    else
    {
        for (j=0; j<ne; j++)
            nr2[j]=0;
        fallos=2*ne;

    }

//delete [] aux;
    //delete [] aux_p;
    //delete [] aux_n;

 if (regla_valida and valida) {

  // Evaluar si es posible eliminar alguna regla ya existente
  int reglas_eliminables=0;
  int reglas_otra_clase=0;
  double aciertos_adicionales=0;
  int ganancia_en_variables=(variables_relevantes+funciones_relevantes);
  int ganancia_funciones_en_variables=funciones_irrelevantes;
  double diferencia_peso=peso;
  int diferencia_estables=variables_estables+funciones_estables;
  int regla_afectada;
  double max_porcen=0;
  double balance;
  //int *lista_reglas_eliminables= new int[nrule];
  int diferencia_ceros;
  double diferencia_n_mas_n_menos;

  n_reglas_a_eliminar=0;
  for (int i=tam_rule_not_removed; i<nrule; i++){
    if ( /*(basica_para_regla[i]==0) or (basica_para_regla[i]-basica_esta_regla[i]==0) or */
       (basica_esta_regla[i]-basica_para_regla[i])+ (errores_para_regla[i]-errores_esta_regla[i])>=0){

      if (true /*reglas_eliminables==0 or R.Class_of_rule(i)!=clase*/) {
        reglas_a_eliminar[n_reglas_a_eliminar]=i;
        n_reglas_a_eliminar++;
        reglas_eliminables++;
        if (R.Class_of_rule(i)!=clase)
          reglas_otra_clase++;

        //if (R.Class_of_rule(i)==clase)
          aciertos_adicionales+=(basica_esta_regla[i]-basica_para_regla[i])+ (errores_para_regla[i]-errores_esta_regla[i]);
        /*else
          aciertos_adicionales+=(errores_para_regla[i]-basica_esta_regla[i])-(basica_para_regla[i]-errores_esta_regla[i]);*/

        ganancia_en_variables=R.Variables_of_rule(i)-ganancia_en_variables;
        //ganancia_en_variables-=R.Variables_of_rule(i);
	diferencia_n_mas_n_menos = R.Get_N_Mas(i) - R.Get_N_Menos(i);
        ganancia_funciones_en_variables=R.New_Atributes_of_rule(i)-ganancia_funciones_en_variables;
        diferencia_peso-=R.Get_Weight(i);
        //diferencia_estables-=R.Stable_of_rule(i);
        diferencia_estables=0;
        balance=Diferencia_Adaptacion_Negativa(nr2,R.SeeAdaptaciones_of_rule(i),ne);
	diferencia_ceros = ceros_variables - R.Numero_Ceros_of_rule(i);
        regla_afectada=i;
      }
    }
    else  /*if (basica_esta_regla[i]>0 and basica_para_regla[i]-basica_esta_regla[i]!=0)*/
             if (R.Class_of_rule(i)==clase and (1.0*basica_esta_regla[i]+errores_para_regla[i])/(basica_para_regla[i]+errores_esta_regla[i])>max_porcen)
               max_porcen= (1.0*basica_esta_regla[i]+errores_para_regla[i])/(basica_para_regla[i]+errores_esta_regla[i]);
  }

  /*if (reglas_eliminables>=2){
   cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n" << endl;
   for (int i=0; i<nrule; i++){
    if ( (basica_esta_regla[i]-basica_para_regla[i])+ (errores_para_regla[i]-errores_esta_regla[i])>=0 ) ){

      cout << "\t\t\t\tRegla: " << i << endl;
      cout << "\t\t\t\t\tbasica esta regla: "<< basica_esta_regla[i] << endl;
      cout << "\t\t\t\t\tbasica para regla: "<< basica_para_regla[i] << endl;
      cout << "\t\t\t\t\terrores para regla: "<< errores_para_regla[i] << endl;
      cout << "\t\t\t\t\terrores esta regla: "<< errores_esta_regla[i] << endl;
      cout << "\t\t\t\t\tbasica(esta-para): "<< basica_esta_regla[i]-basica_para_regla[i] << endl;
      cout << "\t\t\t\t\terrores(para-esta): "<< errores_para_regla[i]-errores_esta_regla[i] << endl;
      cout << "\t\t\t\t\t\t Aciertos Adicionales: "<< (basica_esta_regla[i]-basica_para_regla[i])+ (errores_para_regla[i]-errores_esta_regla[i]) << endl << endl;
      //sleep(2);
    }
   }
   cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n" << endl;
  }*/


  //double granularidad = 2.0/ponderacion_por_clase[clase];

  double granularidad = 0;

  if (ponderacion_por_clase[clase]>0)
    positivos = positivos/ponderacion_por_clase[clase];
  else
    positivos = 0;

  if (aciertos+aciertos_adicionales-fallos>granularidad){
    //cout << "aciertos > fallos (+= " << positivos/ponderacion_por_clase[clase] << ", -= " << negativos << endl;
    if (ponderacion_por_clase[clase]>0) {
      if (positivos-negativos<0)
        valor[0]=-0.9;
      else
        //valor[0]=(positivos/ponderacion_por_clase[clase])-(negativos/(ne-ponderacion_por_clase[clase]));
        //valor[0]=(positivos-negativos)/(ponderacion_por_clase[clase]);
        valor[0]=positivos-negativos;

      if (positivos==0 && menos_infinito){
	  valor[0]=-2*ne;
      }
    }
    else
      valor[0]=-2*ne;
  }
  else if (aciertos+aciertos_adicionales>granularidad){
           if (!reducir_reglas){
             valor[0]=aciertos+aciertos_adicionales-fallos-granularidad;
             if (valor[0]==0)
                valor[0]=-ne;
           }
           else {
             valor[0]=aciertos+aciertos_adicionales-fallos-granularidad;
             if (valor[0]==0)
                valor[0]=-ne;
             else if ( reglas_eliminables>=(-valor[0])+2)
                    valor[0]=aciertos+aciertos_adicionales-fallos;
                 else
                   valor[0]=-2*ne;
           }
       }
       else
           //valor[0]=aciertos-fallos;
           valor[0]=-(ne+fallos);


   valor[6]= aciertos_adicionales;
   valor[5]= aciertos-fallos;
   //valor[4]= variables_estables/*+(0.00001)*unos_variables*/;
   //valor[4]= -ceros_variables;
   valor[4]= 0;
   valor[3]= -(variables_relevantes+funciones_relevantes);
   //valor[2]= -negativos_global;
   valor[2]=afecta;
   //valor[1]= afecta;
   valor[1]= reglas_eliminables;
   /* if (reglas_eliminables>=1 and valor[0]>0){
      valor[2]=valor[0];
      valor[0]=valor[0]+((1.0*aciertos_adicionales)/ponderacion_por_clase[clase]);
      valor[4]=reglas_otra_clase;
   }*/

   if ( (reglas_eliminables>=1 and reglas_eliminables<2) and (valor[0]==-ne or valor[0]==0) ){
     if (!discriminante){
       valor[1]=1;
       //valor[2]=balance;
       valor[2]=(n_mas-n_menos) - diferencia_n_mas_n_menos;
       valor[3]=-ganancia_en_variables;
       //valor[4]=diferencia_ceros;
       valor[4]=0;
       //valor[4]=ganancia_funciones_en_variables;
       //valor[5]=diferencia_estables;
       valor[5]=ceros_variables;
       valor[6]=max_porcen;
     }
     else {
       valor[1]=1;
       //valor[3]=0;
       if (aciertos_adicionales>0)
         valor[2]=0.1+aciertos_adicionales;
       else
         /*if (balance>0 and ganancia_en_variables>0){
           valor[2]=1;
         }
         else {
           valor[2]=0;
         }*/
	 //valor[2]=(n_mas-n_menos) - diferencia_n_mas_n_menos;
	 valor[2]=diferencia_peso;
       valor[3]=-ganancia_en_variables;
       //valor[4]=ganancia_funciones_en_variables;
       valor[4]=0;
       //valor[4]=diferencia_ceros;
       valor[5]=ceros_variables;
       //valor[5]=diferencia_estables;
       //valor[6]=max_porcen;
       valor[6]=-R.Get_Weight(reglas_eliminables);
       //valor[2]=balance;
     }
   }
   //valor[6]= fallos;

   if (esp_gen and positivos_global==0){
    valor[0]=-fallos;
   }
   else if (esp_gen and positivos==0){
     valor[0]=-fallos;
   }
   else if (esp_gen and valor[0]<=-0.9){
     valor[0]=((1.0+(aciertos-fallos))/(aciertos+fallos+1)-1)/2.0-0.25;
   }

   if (esp_gen){
     //valor[2]=-DistanciaHamming(regla, initial_rule);
     valor[1]=aciertos;
     //valor[2]= n_mas;
     valor[2]=0;
     //valor[3]= -n_menos;
   }

   //if (reglas_eliminables>=2)
     //valor[3]=-valor[3];

   //delete [] lista_reglas_eliminables;
  }
  else {
    valor[0]=valor[1]=valor[2]=valor[3]=valor[4]=valor[5]=valor[6]=-2*ne;
    peso=0;
  }


  if (nrule>0 and activada_memoria){
    delete [] basica_esta_regla;
    delete [] errores_esta_regla;
  }

  delete [] ejemplos_afectado_por_esta_regla;
  delete [] nr2;

  return peso;
}



double Bondad_E (multipoblacion &G, const VectorVar &V, const example_set &E, ruleset &R, double *valor, int N_valor,
               int elemento_i, bool discriminante, bool reducir_reglas, int *ponderacion_por_clase, int nrule,
               double ***Tabla, double ***Tabla_Nuevas, VectorVar &F, int *basica_para_regla, int *errores_para_regla, int *reglas_a_eliminar, int &n_reglas_a_eliminar,
	       bool esp_gen, string initial_rule, double nmas, double nmenos)
{

    int b1, b2, n1, n2, r1, r2, r3;
    char *nb1, *nb2;
    int *nn1,*nn2,r;
    double *nr1, *nr2, *nr3;
    string regla, regla_funciones;
    int var1, var2;
    int ne=E.N_Examples();
    vectordouble w(V.N_Antecedente());
    double positivos=0, negativos=0, rango, eficacia, negativos_que_conllevan_error=0, positivos_a_incrementar=0,
                                  positivos_sobrecubrimiento=0, negativos2=0;
    int j;
    double aciertos=0, fallos=0, aciertos_aqui=0, fallos_aqui=0, sin_cubrir=0, aciertos_no_cubiertos=0, max_pos=0,
                              max_neg=0, max_neg2=0;
    int clase = G.ClaseIndividuo(elemento_i);
    bool menos_infinito=true, valida=true;
    int variables_irrelevantes=0, funciones_irrelevantes=0, variables_estables=0, funciones_estables=0,
        ceros_variables=0, ceros_funciones=0, unos_variables=0, unos_funciones=0,
	variables_relevantes, funciones_relevantes;

    nrule= R.N_rule();


    nb1=G.Subpoblacion_Binaria(0,elemento_i,b1);
    char *s= new char[b1+1];
    valor[2]=0;  //numero de unos en el cromosoma
    for (j=0; j<b1; j++)
    {
        s[j]=nb1[j];
        if (nb1[j]=='1')
            valor[2]++;
        s[b1]='\0';
    }
    regla=s;
    delete []s;


    nb2=G.Subpoblacion_Binaria(1,elemento_i,b2);
    s= new char[b2+1];
    for (j=0; j<b2; j++)
    {
        s[j]=nb2[j];
        if (nb2[j]=='1')
            valor[2]++;
        s[b2]='\0';
    }
    regla_funciones=s;
    delete []s;


    nn1=G.Subpoblacion_Entera(0,elemento_i,n1);
    nn2=G.Subpoblacion_Entera(1,elemento_i,n2);
    nr1=G.Subpoblacion_Real(0,elemento_i,r1);
    nr2 = new double[ne];


    valor[1]=0; //numero de variable irrelevantes
    valor[2]=0; //funciones irrelevantes
    valor[3]=0; //asignaciones comprensibles
    valor[4]=0; //numero de ceros
    int pos=0; // posicion en el cromosoma nivel de valor
    int tamanio; // numero de genes de una variable

    //cout << "Umbral: " << nr1[r1-1] << endl;
    for (j=0; j<r1-1; j++)
    {
        tamanio= V.Variable(j).N_etiquetas();
        if (nr1[j]<nr1[r1-1])
            variables_irrelevantes++;
        else
        {
           if (Relevante(&(nb1[pos]),tamanio)){
              if (Asignacion_Comprensible(&(nb1[pos]),tamanio)){
                variables_estables++;
                //ceros_variables+=NumeroCeros(&(nb1[pos]),tamanio);
                //unos_variables+=NumeroUnos(&(nb1[pos]),tamanio);
              }
              //else {
                ceros_variables+=NumeroCeros(&(nb1[pos]),tamanio);
                unos_variables+=NumeroUnos(&(nb1[pos]),tamanio);
                // Asignacion no comprensible
                // considero entonces la regla como invalida
                //valida=false;
              //}
           }
           else {
              nr1[j]=nr1[r1-1]-0.001;
              variables_irrelevantes++;
           }

        }
        pos=pos+tamanio;
    }

    // Quitar funciones repetidas

    for (int i=1; i<n2; i++){
      if (nn2[i]!=0) {
        int t=0;
        while (t<i and nn2[t]!=nn2[i])
          t++;
        if (t<i)
          nn2[i]=0;
      }
    }

    variables_relevantes = (r1-1) - variables_irrelevantes;


    pos=0;
    for (j=0; j<n2; j++)
    {
        if (nn2[j]!=0)
        {
           tamanio= F.SizeDomain(nn2[j]-1);
           if (Relevante(&(nb2[pos]),tamanio)){
              if (Asignacion_Comprensible(&(nb2[pos]),tamanio)){
                 funciones_estables++;
                 ceros_funciones+=NumeroCeros(&(nb2[pos]),tamanio);
                 unos_funciones+=NumeroUnos(&(nb2[pos]),tamanio);
              }
              else {
                // Asignacion no comprensible
                // considero entonces la regla como invalida
                //valida=false;
              }
           }
           else {
              funciones_irrelevantes++;
              nn2[j]=0;
           }
        }
        else {
            tamanio=2;
            funciones_irrelevantes++;
        }
        tamanio=2;
        pos=pos+tamanio;
    }

    funciones_relevantes = (n2) - funciones_irrelevantes;


//Distingue entre descripciones discriminantes o descriptivas
    /*if (!discriminante)
        valor[1]=r1-1+n2-valor[1]-valor[2];*/


    bool activada_memoria=false;
    int *basica_esta_regla = 0;
    int *errores_esta_regla = 0;
    if (nrule>0) {
      basica_esta_regla = new int[nrule];
      errores_esta_regla = new int[nrule];
      activada_memoria=true;
      for (int i=0; i<nrule; i++){
        basica_esta_regla[i]=0;
        errores_esta_regla[i]=0;
      }
    }

    int *ejemplos_afectado_por_esta_regla = new int[ne];

//double *aux = new double[ne];
    double aux_p;
    double aux_n;
    double peso, kk;
    double fitness0_global_regla=-ne;
// Comprobar que es una regla valida
    bool esta_cubierto;
//bool regla_valida = V.Es_Valida(regla,nr1,nr1[r1-1],valor[2]);
    bool regla_valida = V.Es_Valida(regla,nr1,nr1[r1-1],kk);
    double afecta=0;


    if (valida and !esp_gen)
       valida = ponderacion_por_clase[clase]!=0;
    else {
       valida = true;
       if (!regla_valida){
	 //cout << "Regla considerada como no valida\n";
       }
    }

// Calculo previo para obtener el peso de la regla sobre atributos originales
    if (regla_valida and valida)
    {
        for (j=0; j<ne; j++)
        {

            T_AdaptacionC(V,Tabla,j,nn1[0],aux_p,aux_n);

            //esta_cubierto = E.Is_Covered(j);

            nr2[j]=T_Adaptacion_Antecedente(V,Tabla,j,regla,nr1,nr1[r1-1]);
            //nr2[j]=1;
            nr2[j]=T_Adaptacion_Antecedente_Nuevas(F, Tabla_Nuevas, j, regla_funciones, nn2, n2, nr2[j]);

            if (nr2[j]*aux_p>0)
              afecta++;

            positivos=positivos+(nr2[j]*aux_p);
            negativos=negativos+(nr2[j]*aux_n);
            if (aux_p<aux_n)
                nr2[j]=-nr2[j];

        }

        //for (int i=0; i<nrule; i++)
          //cout <<"Regla " << i <<": " << basica_para_regla[i] << endl;
        //fitness0_global_regla = (positivos_global-negativos_global)/(ponderacion_por_clase[clase]);
    }


//peso = (positivos+1) / (positivos+negativos+1);
//peso = 1.0;
    nmas+=positivos;
    nmenos+=negativos;
    peso = (nmas+1) / (nmas+nmenos+1);
//peso = nr3[0];

 if (positivos>0)
   menos_infinito=false;


 /*if (afecta< 0.20*ponderacion_por_clase[clase])
   valida=false;*/

 if (ponderacion_por_clase[clase]>0)
   afecta = afecta / ponderacion_por_clase[clase];

 double positivos_global=0, negativos_global=0;
// Dado el peso de la regla obtener la bondad de la misma.
    int positivos_aux=0;
    if (regla_valida and valida)
    {
        positivos=0;
        negativos=0;
        for (j=0; j<ne; j++)
        {
            if (nr2[j]>=0)
              positivos_global+= (nr2[j]*peso);
            else
              negativos_global+= (-nr2[j]*peso);

            ejemplos_afectado_por_esta_regla[j]=0;

            esta_cubierto = E.Is_Covered(j);
            if (!esta_cubierto)
            {
                if (nr2[j]>0)
                {
                    if ( (nr2[j]*peso>E.Grade_Is_Negative_Covered(j) ) or (nr2[j]*peso==E.Grade_Is_Negative_Covered(j) and peso>E.Peso_Is_Negative_Covered(j) ) )
                    {
                        positivos=positivos+(nr2[j]*peso);
                        aciertos++;
                        ejemplos_afectado_por_esta_regla[j]=1;
                        if (nrule>0 and E.Rule_Aislada(j)==-1 and E.Rule_Used(j)>=0 and E.Rule_Used(j)<nrule)
                          errores_esta_regla[E.Rule_Used(j)]++;
                    }
                    else {
                        positivos_a_incrementar += (nr2[j]*peso);
                    }
                }
                else if (nr2[j]<0) {

                       if (  (-nr2[j]*peso>E.Grade_Is_Positive_Covered(j)) or (-nr2[j]*peso==E.Grade_Is_Positive_Covered(j) and peso>E.Peso_Is_Positive_Covered(j)) ){
                          //negativos=negativos+(-nr2[j]*peso);
                          //if (E.Grade_Is_Positive_Covered(j)>0)
                            negativos=negativos+((-nr2[j]*peso)/ponderacion_por_clase[E.Class_of_example(V.Consecuente(),j)]);
                          //negativos=negativos+((-nr2[j]*peso)/ponderacion_por_clase[E.Class_of_example(V.Consecuente(),j)]);
                          ejemplos_afectado_por_esta_regla[j]=-1;

                          if (nrule>0 and E.Rule_Aislada(j)==-1 and E.Rule_Used(j)>=0 and E.Rule_Used(j)<nrule)
                            errores_esta_regla[E.Rule_Used(j)]++;
                       }

                }
            }
            else
            {
                if (nr2[j]>0)
                {
                   if (nr2[j]*peso>E.Grade_Is_Negative_Covered(j) or
                      (nr2[j]*peso==E.Grade_Is_Negative_Covered(j) and peso>E.Peso_Is_Negative_Covered(j) ) ){
                        ejemplos_afectado_por_esta_regla[j]=1;
                        if (nrule>0 and E.Rule_Aislada(j)==1 and E.Rule_Used(j)>=0 and E.Rule_Used(j)<nrule)
                           basica_esta_regla[E.Rule_Used(j)]++;
                        positivos_sobrecubrimiento+=1;
                   }

                }
                else  if (nr2[j]<0 and ( (-nr2[j])*peso>E.Grade_Is_Positive_Covered(j) or
                         ((-nr2[j])*peso==E.Grade_Is_Positive_Covered(j) and peso > E.Peso_Is_Positive_Covered(j) ) ) )
                {
                    ejemplos_afectado_por_esta_regla[j]=-1;
                    //negativos=negativos+(-nr2[j]*peso);
                    negativos=negativos+((-nr2[j]*peso)/ponderacion_por_clase[E.Class_of_example(V.Consecuente(),j)]);
                    fallos++;
                    if (nrule>0 and E.Rule_Aislada(j)==1 and E.Rule_Used(j)>=0 and E.Rule_Used(j)<nrule)
                      basica_esta_regla[E.Rule_Used(j)]++;
                }
            }
        }
    }
// La regla es invalida
    else
    {
        for (j=0; j<ne; j++)
            nr2[j]=0;
        fallos=2*ne;

    }

//delete [] aux;
    //delete [] aux_p;
    //delete [] aux_n;

 if (regla_valida and valida) {

  // Evaluar si es posible eliminar alguna regla ya existente
  int reglas_eliminables=0;
  int reglas_otra_clase=0;
  double aciertos_adicionales=0;
  int ganancia_en_variables=(variables_relevantes+funciones_relevantes);
  int ganancia_funciones_en_variables=funciones_irrelevantes;
  double diferencia_peso=peso;
  int diferencia_estables=variables_estables+funciones_estables;
  int regla_afectada;
  double max_porcen=0;
  double balance;
  //int *lista_reglas_eliminables= new int[nrule];
  int diferencia_ceros;
  double diferencia_n_mas_n_menos;

  n_reglas_a_eliminar=0;
  for (int i=0; i<nrule; i++){
    if ( /*(basica_para_regla[i]==0) or (basica_para_regla[i]-basica_esta_regla[i]==0) or */
       (basica_esta_regla[i]-basica_para_regla[i])+ (errores_para_regla[i]-errores_esta_regla[i])>=0){

      if (true /*reglas_eliminables==0 or R.Class_of_rule(i)!=clase*/) {
        reglas_a_eliminar[n_reglas_a_eliminar]=i;
        n_reglas_a_eliminar++;
        reglas_eliminables++;
        if (R.Class_of_rule(i)!=clase)
          reglas_otra_clase++;

        //if (R.Class_of_rule(i)==clase)
          aciertos_adicionales+=(basica_esta_regla[i]-basica_para_regla[i])+ (errores_para_regla[i]-errores_esta_regla[i]);
        /*else
          aciertos_adicionales+=(errores_para_regla[i]-basica_esta_regla[i])-(basica_para_regla[i]-errores_esta_regla[i]);*/

        ganancia_en_variables=R.Variables_of_rule(i)-ganancia_en_variables;
        //ganancia_en_variables-=R.Variables_of_rule(i);
	diferencia_n_mas_n_menos = R.Get_N_Mas(i) - R.Get_N_Menos(i);
        ganancia_funciones_en_variables=R.New_Atributes_of_rule(i)-ganancia_funciones_en_variables;
        diferencia_peso-=R.Get_Weight(i);
        //diferencia_estables-=R.Stable_of_rule(i);
        diferencia_estables=0;
        balance=Diferencia_Adaptacion_Negativa(nr2,R.SeeAdaptaciones_of_rule(i),ne);
	diferencia_ceros = ceros_variables - R.Numero_Ceros_of_rule(i);
        regla_afectada=i;
      }
    }
    else  /*if (basica_esta_regla[i]>0 and basica_para_regla[i]-basica_esta_regla[i]!=0)*/
             if (R.Class_of_rule(i)==clase and (1.0*basica_esta_regla[i]+errores_para_regla[i])/(basica_para_regla[i]+errores_esta_regla[i])>max_porcen)
               max_porcen= (1.0*basica_esta_regla[i]+errores_para_regla[i])/(basica_para_regla[i]+errores_esta_regla[i]);
  }

  /*if (reglas_eliminables>=2){
   cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n" << endl;
   for (int i=0; i<nrule; i++){
    if ( (basica_esta_regla[i]-basica_para_regla[i])+ (errores_para_regla[i]-errores_esta_regla[i])>=0 ) ){

      cout << "\t\t\t\tRegla: " << i << endl;
      cout << "\t\t\t\t\tbasica esta regla: "<< basica_esta_regla[i] << endl;
      cout << "\t\t\t\t\tbasica para regla: "<< basica_para_regla[i] << endl;
      cout << "\t\t\t\t\terrores para regla: "<< errores_para_regla[i] << endl;
      cout << "\t\t\t\t\terrores esta regla: "<< errores_esta_regla[i] << endl;
      cout << "\t\t\t\t\tbasica(esta-para): "<< basica_esta_regla[i]-basica_para_regla[i] << endl;
      cout << "\t\t\t\t\terrores(para-esta): "<< errores_para_regla[i]-errores_esta_regla[i] << endl;
      cout << "\t\t\t\t\t\t Aciertos Adicionales: "<< (basica_esta_regla[i]-basica_para_regla[i])+ (errores_para_regla[i]-errores_esta_regla[i]) << endl << endl;
      //sleep(2);
    }
   }
   cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n" << endl;
  }*/


  //double granularidad = 2.0/ponderacion_por_clase[clase];

  double granularidad = 0;

  if (ponderacion_por_clase[clase]>0)
    positivos = positivos/ponderacion_por_clase[clase];
  else
    positivos = 0;

  if (aciertos+aciertos_adicionales-fallos>granularidad){
    //cout << "aciertos > fallos (+= " << positivos/ponderacion_por_clase[clase] << ", -= " << negativos << endl;
    if (ponderacion_por_clase[clase]>0) {
      if (positivos-negativos<0)
        valor[0]=-0.9;
      else
        //valor[0]=(positivos/ponderacion_por_clase[clase])-(negativos/(ne-ponderacion_por_clase[clase]));
        //valor[0]=(positivos-negativos)/(ponderacion_por_clase[clase]);
        valor[0]=positivos-negativos;

      if (positivos==0 && menos_infinito){
	  valor[0]=-2*ne;
      }
    }
    else
      valor[0]=-2*ne;
  }
  else if (aciertos+aciertos_adicionales>granularidad){
           if (!reducir_reglas){
             valor[0]=aciertos+aciertos_adicionales-fallos-granularidad;
             if (valor[0]==0)
                valor[0]=-ne;
           }
           else {
             valor[0]=aciertos+aciertos_adicionales-fallos-granularidad;
             if (valor[0]==0)
                valor[0]=-ne;
             else if ( reglas_eliminables>=(-valor[0])+2)
                    valor[0]=aciertos+aciertos_adicionales-fallos;
                 else
                   valor[0]=-2*ne;
           }
       }
       else
           //valor[0]=aciertos-fallos;
           valor[0]=-(ne+fallos);


   valor[6]= aciertos_adicionales;
   valor[5]= aciertos-fallos;
   //valor[4]= variables_estables/*+(0.00001)*unos_variables*/;
   //valor[4]= -ceros_variables;
   valor[4]= 0;
   valor[3]= -(variables_relevantes+funciones_relevantes);
   //valor[2]= -negativos_global;
   valor[2]=afecta;
   //valor[1]= afecta;
   valor[1]= reglas_eliminables;
   /* if (reglas_eliminables>=1 and valor[0]>0){
      valor[2]=valor[0];
      valor[0]=valor[0]+((1.0*aciertos_adicionales)/ponderacion_por_clase[clase]);
      valor[4]=reglas_otra_clase;
   }*/

   if ( (reglas_eliminables>=1 and reglas_eliminables<2) and (valor[0]==-ne or valor[0]==0) ){
     if (!discriminante){
       valor[1]=1;
       //valor[2]=balance;
       valor[2]=(nmas-nmenos) - diferencia_n_mas_n_menos;
       valor[3]=-ganancia_en_variables;
       //valor[4]=diferencia_ceros;
       valor[4]=0;
       //valor[4]=ganancia_funciones_en_variables;
       //valor[5]=diferencia_estables;
       valor[5]=ceros_variables;
       valor[6]=max_porcen;
     }
     else {
       valor[1]=1;
       //valor[3]=0;
       if (aciertos_adicionales>0)
         valor[2]=0.1+aciertos_adicionales;
       else
         /*if (balance>0 and ganancia_en_variables>0){
           valor[2]=1;
         }
         else {
           valor[2]=0;
         }*/
	 //valor[2]=(n_mas-n_menos) - diferencia_n_mas_n_menos;
	 valor[2]=diferencia_peso;
       valor[3]=-ganancia_en_variables;
       //valor[4]=ganancia_funciones_en_variables;
       valor[4]=0;
       //valor[4]=diferencia_ceros;
       valor[5]=ceros_variables;
       //valor[5]=diferencia_estables;
       //valor[6]=max_porcen;
       valor[6]=-R.Get_Weight(reglas_eliminables);
       //valor[2]=balance;
     }
   }
   //valor[6]= fallos;

   if (esp_gen and positivos_global==0){
    valor[0]=-fallos;
   }
   else if (esp_gen and positivos==0){
     valor[0]=-fallos;
   }
   else if (esp_gen and valor[0]<=-0.9){
     valor[0]=((1.0+(aciertos-fallos))/(aciertos+fallos+1)-1)/2.0-0.25;
   }

   if (esp_gen){
     //valor[2]=-DistanciaHamming(regla, initial_rule);
     valor[1]=aciertos;
     //valor[2]= n_mas;
     valor[2]=0;
     //valor[3]= -n_menos;
   }

   //if (reglas_eliminables>=2)
     //valor[3]=-valor[3];

   //delete [] lista_reglas_eliminables;
  }
  else {
    valor[0]=valor[1]=valor[2]=valor[3]=valor[4]=valor[5]=valor[6]=-2*ne;
    peso=0;
  }


  if (nrule>0 and activada_memoria){
    delete [] basica_esta_regla;
    delete [] errores_esta_regla;
  }

  delete [] ejemplos_afectado_por_esta_regla;
  delete [] nr2;

  return peso;
}




double Bondad_Final (multipoblacion &G, const VectorVar &V, const example_set &E, ruleset &R, double *valor, int N_valor,
               int elemento_i, bool discriminante, bool reducir_reglas, int *ponderacion_por_clase, int nrule,
               double ***Tabla, double ***Tabla_Nuevas, VectorVar &F, int *basica_para_regla, int *errores_para_regla, int *reglas_a_eliminar, int &n_reglas_a_eliminar)
{

    int b1, b2, n1, n2, r1, r2, r3;
    char *nb1, *nb2;
    int *nn1,*nn2,r;
    double *nr1, *nr2, *nr3;
    string regla, regla_funciones;
    int var1, var2;
    int ne=E.N_Examples();
    vectordouble w(V.N_Antecedente());
    double positivos=0, negativos=0, rango, eficacia, negativos_que_conllevan_error=0, positivos_a_incrementar=0,
                                  positivos_sobrecubrimiento=0, negativos2=0;
    int j;
    double aciertos=0, fallos=0, aciertos_aqui=0, fallos_aqui=0, sin_cubrir=0, aciertos_no_cubiertos=0, max_pos=0,
                              max_neg=0, max_neg2=0;
    int clase = G.ClaseIndividuo(elemento_i);
    bool menos_infinito=true, valida=true;
    int variables_irrelevantes=0, funciones_irrelevantes=0, variables_estables=0, funciones_estables=0, ceros_variables=0, ceros_funciones=0, unos_variables=0, unos_funciones=0;


    nb1=G.Subpoblacion_Binaria(0,elemento_i,b1);
    char *s= new char[b1+1];
    valor[2]=0;  //numero de unos en el cromosoma
    for (j=0; j<b1; j++)
    {
        s[j]=nb1[j];
        if (nb1[j]=='1')
            valor[2]++;
        s[b1]='\0';
    }
    regla=s;
    delete []s;


    nb2=G.Subpoblacion_Binaria(1,elemento_i,b2);
    s= new char[b2+1];
    for (j=0; j<b2; j++)
    {
        s[j]=nb2[j];
        if (nb2[j]=='1')
            valor[2]++;
        s[b2]='\0';
    }
    regla_funciones=s;
    delete []s;


    nn1=G.Subpoblacion_Entera(0,elemento_i,n1);
    nn2=G.Subpoblacion_Entera(1,elemento_i,n2);
    nr1=G.Subpoblacion_Real(0,elemento_i,r1);
    nr2 = new double[ne];


    valor[1]=0; //numero de variable irrelevantes
    valor[2]=0; //funciones irrelevantes
    valor[3]=0; //asignaciones comprensibles
    valor[4]=0; //numero de ceros
    int pos=0; // posicion en el cromosoma nivel de valor
    int tamanio; // numero de genes de una variable
    for (j=0; j<r1-1; j++)
    {
        tamanio= V.Variable(j).N_etiquetas();
        if (nr1[j]<nr1[r1-1])
            variables_irrelevantes++;
        else
        {
           if (Relevante(&(nb1[pos]),tamanio)){
              if (Asignacion_Comprensible(&(nb1[pos]),tamanio)){
                variables_estables++;
                ceros_variables+=NumeroCeros(&(nb1[pos]),tamanio);
                unos_variables+=NumeroUnos(&(nb1[pos]),tamanio);
              }
              else {
                // Asignacion no comprensible
                // considero entonces la regla como invalida
                //valida=false;
              }
           }
           else {
              nr1[j]=nr1[r1-1]-0.001;
              variables_irrelevantes++;
           }

        }
        pos=pos+tamanio;
    }

    // Quitar funciones repetidas

    for (int i=1; i<n2; i++){
      if (nn2[i]!=0) {
        int t=0;
        while (t<i and nn2[t]!=nn2[i])
          t++;
        if (t<i)
          nn2[i]=0;
      }
    }


    pos=0;
    for (j=0; j<n2; j++)
    {
        if (nn2[j]!=0)
        {
           tamanio= F.SizeDomain(nn2[j]-1);
           if (Relevante(&(nb2[pos]),tamanio)){
              if (Asignacion_Comprensible(&(nb2[pos]),tamanio)){
                 funciones_estables++;
                 ceros_funciones+=NumeroCeros(&(nb2[pos]),tamanio);
                 unos_funciones+=NumeroUnos(&(nb2[pos]),tamanio);
              }
              else {
                // Asignacion no comprensible
                // considero entonces la regla como invalida
                //valida=false;
              }
           }
           else {
              funciones_irrelevantes++;
              nn2[j]=0;
           }
        }
        else {
            tamanio=2;
            funciones_irrelevantes++;
        }
        tamanio=2;
        pos=pos+tamanio;
    }



//Distingue entre descripciones discriminantes o descriptivas
    /*if (!discriminante)
        valor[1]=r1-1+n2-valor[1]-valor[2];*/



    int *basica_esta_regla = 0;
    int *errores_esta_regla = 0;
    if (nrule>0) {
      basica_esta_regla = new int[nrule];
      errores_esta_regla = new int[nrule];
      for (int i=0; i<nrule; i++){
        basica_esta_regla[i]=0;
        errores_esta_regla[i]=0;
      }
    }

    int *ejemplos_afectado_por_esta_regla = new int[ne];

//double *aux = new double[ne];
    double aux_p;
    double aux_n;
    double peso, kk;
    double fitness0_global_regla=-ne;
// Comprobar que es una regla valida
    bool esta_cubierto;
//bool regla_valida = V.Es_Valida(regla,nr1,nr1[r1-1],valor[2]);
    bool regla_valida = V.Es_Valida(regla,nr1,nr1[r1-1],kk);
    double afecta=0;

    if (valida)
       valida = ponderacion_por_clase[clase]!=0;

// Calculo previo para obtener el peso de la regla sobre atributos originales
    if (regla_valida and valida)
    {
        for (j=0; j<ne; j++)
        {

            T_AdaptacionC(V,Tabla,j,nn1[0],aux_p,aux_n);

            //esta_cubierto = E.Is_Covered(j);

            nr2[j]=T_Adaptacion_Antecedente(V,Tabla,j,regla,nr1,nr1[r1-1]);
            //nr2[j]=1;
            nr2[j]=T_Adaptacion_Antecedente_Nuevas(F, Tabla_Nuevas, j, regla_funciones, nn2, n2, nr2[j]);

            if (nr2[j]*aux_p>0)
              afecta++;

            positivos=positivos+(nr2[j]*aux_p);
            negativos=negativos+(nr2[j]*aux_n);
            if (aux_p<aux_n)
                nr2[j]=-nr2[j];

        }

        //for (int i=0; i<nrule; i++)
          //cout <<"Regla " << i <<": " << basica_para_regla[i] << endl;
        //fitness0_global_regla = (positivos_global-negativos_global)/(ponderacion_por_clase[clase]);
    }


//peso = (positivos+1) / (positivos+negativos+1);
//peso = 1.0;
    peso = (positivos+1) / (positivos+negativos+1);
//peso = nr3[0];

 if (positivos>0)
   menos_infinito=false;


 /*if (afecta< 0.20*ponderacion_por_clase[clase])
   valida=false;*/

 if (ponderacion_por_clase[clase]>0)
   afecta = afecta / ponderacion_por_clase[clase];

 double positivos_global=0, negativos_global=0;
// Dado el peso de la regla obtener la bondad de la misma.
    int positivos_aux=0;
    if (regla_valida and valida)
    {
        positivos=0;
        negativos=0;
        for (j=0; j<ne; j++)
        {
            if (nr2[j]>=0)
              positivos_global+= (nr2[j]*peso);
            else
              negativos_global+= (-nr2[j]*peso);

            ejemplos_afectado_por_esta_regla[j]=0;

            esta_cubierto = E.Is_Covered(j);
            if (!esta_cubierto)
            {
                if (nr2[j]>0)
                {
                    if ( (nr2[j]*peso>E.Grade_Is_Negative_Covered(j) ) or (nr2[j]*peso==E.Grade_Is_Negative_Covered(j) and peso>E.Peso_Is_Negative_Covered(j) ) )
                    {
                        positivos=positivos+(nr2[j]*peso);
                        aciertos++;
                        ejemplos_afectado_por_esta_regla[j]=1;
                        if (nrule>0 and E.Rule_Aislada(j)==-1 and E.Rule_Used(j)>=0)
                          errores_esta_regla[E.Rule_Used(j)]++;
                    }
                    else {
                        positivos_a_incrementar += (nr2[j]*peso);
                    }
                }
                else if (nr2[j]<0) {

                       if (  (-nr2[j]*peso>E.Grade_Is_Positive_Covered(j)) or (-nr2[j]*peso==E.Grade_Is_Positive_Covered(j) and peso>E.Peso_Is_Positive_Covered(j)) ){
                          //negativos=negativos+(-nr2[j]*peso);
                          //if (E.Grade_Is_Positive_Covered(j)>0)
                            negativos=negativos+((-nr2[j]*peso)/ponderacion_por_clase[E.Class_of_example(V.Consecuente(),j)]);
                          //negativos=negativos+((-nr2[j]*peso)/ponderacion_por_clase[E.Class_of_example(V.Consecuente(),j)]);
                          ejemplos_afectado_por_esta_regla[j]=-1;

                          if (nrule>0 and E.Rule_Aislada(j)==-1 and E.Rule_Used(j)>=0)
                            errores_esta_regla[E.Rule_Used(j)]++;
                       }

                }
            }
            else
            {
                if (nr2[j]>0)
                {
                   if (nr2[j]*peso>E.Grade_Is_Negative_Covered(j) or
                      (nr2[j]*peso==E.Grade_Is_Negative_Covered(j) and peso>E.Peso_Is_Negative_Covered(j) ) ){
                        ejemplos_afectado_por_esta_regla[j]=1;
                        if (nrule>0 and E.Rule_Aislada(j)==1 and E.Rule_Used(j)>=0)
                           basica_esta_regla[E.Rule_Used(j)]++;
                        positivos_sobrecubrimiento+=1;
                   }

                }
                else  if (nr2[j]<0 and ( (-nr2[j])*peso>E.Grade_Is_Positive_Covered(j) or
                         ((-nr2[j])*peso==E.Grade_Is_Positive_Covered(j) and peso > E.Peso_Is_Positive_Covered(j) ) ) )
                {
                    ejemplos_afectado_por_esta_regla[j]=-1;
                    //negativos=negativos+(-nr2[j]*peso);
                    negativos=negativos+((-nr2[j]*peso)/ponderacion_por_clase[E.Class_of_example(V.Consecuente(),j)]);
                    fallos++;
                    if (nrule>0 and E.Rule_Aislada(j)==1 and E.Rule_Used(j)>=0)
                      basica_esta_regla[E.Rule_Used(j)]++;
                }
            }
        }
    }
// La regla es invalida
    else
    {
        for (j=0; j<ne; j++)
            nr2[j]=0;
        fallos=2*ne;

    }

//delete [] aux;
    //delete [] aux_p;
    //delete [] aux_n;

 if (regla_valida and valida) {

  // Evaluar si es posible eliminar alguna regla ya existente
  int reglas_eliminables=0;
  int reglas_otra_clase=0;
  double aciertos_adicionales=0;
  int ganancia_en_variables=(variables_irrelevantes+funciones_irrelevantes);
  int ganancia_funciones_en_variables=funciones_irrelevantes;
  double diferencia_peso=peso;
  int diferencia_estables=variables_estables+funciones_estables;
  int regla_afectada;
  double max_porcen=0;
  double balance;
  int *lista_reglas_eliminables= new int[nrule];

  n_reglas_a_eliminar=0;
  for (int i=0; i<nrule; i++){
    if ( /*(basica_para_regla[i]==0) or (basica_para_regla[i]-basica_esta_regla[i]==0) or */
       (basica_esta_regla[i]-basica_para_regla[i])+ (errores_para_regla[i]-errores_esta_regla[i])>=0){

      if (true /*reglas_eliminables==0 or R.Class_of_rule(i)!=clase*/) {
        reglas_a_eliminar[n_reglas_a_eliminar]=i;
        n_reglas_a_eliminar++;
        reglas_eliminables++;
        if (R.Class_of_rule(i)!=clase)
          reglas_otra_clase++;

        //if (R.Class_of_rule(i)==clase)
          aciertos_adicionales+=(basica_esta_regla[i]-basica_para_regla[i])+ (errores_para_regla[i]-errores_esta_regla[i]);
        /*else
          aciertos_adicionales+=(errores_para_regla[i]-basica_esta_regla[i])-(basica_para_regla[i]-errores_esta_regla[i]);*/

        ganancia_en_variables-=R.Variables_of_rule(i);
        ganancia_funciones_en_variables-=R.New_Atributes_of_rule(i);
        diferencia_peso-=R.Get_Weight(i);
        diferencia_estables-=R.Stable_of_rule(i);
        balance=Diferencia_Adaptacion_Negativa(nr2,R.SeeAdaptaciones_of_rule(i),ne);
        regla_afectada=i;
      }
      else {
        lista_reglas_eliminables[reglas_eliminables]=i;
        reglas_eliminables++;
        int incremento_temporal_aciertos= aciertos_adicionales;
        aciertos_adicionales=(basica_esta_regla[i]-basica_para_regla[i])+ (errores_para_regla[i]-errores_esta_regla[i]);
        double *adapt_nueva_eliminada;
        int j=0;
        while (j<ne and aciertos_adicionales>=0){
          if (ejemplos_afectado_por_esta_regla[j]==0 and E.Rule_Aislada(j)>0 and E.Rule_Aislada(j)<=reglas_eliminables-reglas_otra_clase){
            int k=0, tope=E.Rule_Aislada(j), cont=0;
            while (k<reglas_eliminables and cont+(reglas_eliminables-reglas_otra_clase-k)>=tope){
              if (R.Class_of_rule(lista_reglas_eliminables[k])==clase){
                adapt_nueva_eliminada = R.SeeAdaptaciones_of_rule(lista_reglas_eliminables[k]);
                if (adapt_nueva_eliminada[j]>E.Grade_Is_Negative_Covered(j) or (adapt_nueva_eliminada[j]>0 and adapt_nueva_eliminada[j]==E.Grade_Is_Negative_Covered(j) and R.Get_Weight(lista_reglas_eliminables[k])>=E.Peso_Is_Negative_Covered(j))    ) {
                  cont++;
                }
              }
              k++;
            }
            if (cont==tope)
                aciertos_adicionales--;
          }
          j++;
        }

        if (aciertos_adicionales<0){
          reglas_eliminables--;
          aciertos_adicionales=incremento_temporal_aciertos;
        }
        else {
          aciertos_adicionales+=incremento_temporal_aciertos;
        }

      }

    }
    else  /*if (basica_esta_regla[i]>0 and basica_para_regla[i]-basica_esta_regla[i]!=0)*/
             if (R.Class_of_rule(i)==clase and (1.0*basica_esta_regla[i]+errores_para_regla[i])/(basica_para_regla[i]+errores_esta_regla[i])>max_porcen)
               max_porcen= (1.0*basica_esta_regla[i]+errores_para_regla[i])/(basica_para_regla[i]+errores_esta_regla[i]);
  }

  /*if (reglas_eliminables>=2){
   cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n" << endl;
   for (int i=0; i<nrule; i++){
    if ( (basica_esta_regla[i]-basica_para_regla[i])+ (errores_para_regla[i]-errores_esta_regla[i])>=0 ) ){

      cout << "\t\t\t\tRegla: " << i << endl;
      cout << "\t\t\t\t\tbasica esta regla: "<< basica_esta_regla[i] << endl;
      cout << "\t\t\t\t\tbasica para regla: "<< basica_para_regla[i] << endl;
      cout << "\t\t\t\t\terrores para regla: "<< errores_para_regla[i] << endl;
      cout << "\t\t\t\t\terrores esta regla: "<< errores_esta_regla[i] << endl;
      cout << "\t\t\t\t\tbasica(esta-para): "<< basica_esta_regla[i]-basica_para_regla[i] << endl;
      cout << "\t\t\t\t\terrores(para-esta): "<< errores_para_regla[i]-errores_esta_regla[i] << endl;
      cout << "\t\t\t\t\t\t Aciertos Adicionales: "<< (basica_esta_regla[i]-basica_para_regla[i])+ (errores_para_regla[i]-errores_esta_regla[i]) << endl << endl;
      //sleep(2);
    }
   }
   cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n" << endl;
  }*/


  //double granularidad = 2.0/ponderacion_por_clase[clase];

  double granularidad = 0;

  if (ponderacion_por_clase[clase]>0)
    positivos = positivos/ponderacion_por_clase[clase];

  if (aciertos+aciertos_adicionales-fallos>granularidad){
    //cout << "aciertos > fallos (+= " << positivos/ponderacion_por_clase[clase] << ", -= " << negativos << endl;
    if (ponderacion_por_clase[clase]>0) {
      if (positivos-negativos<0)
        valor[0]=-0.9;
      else
        //valor[0]=(positivos/ponderacion_por_clase[clase])-(negativos/(ne-ponderacion_por_clase[clase]));
        //valor[0]=(positivos-negativos)/(ponderacion_por_clase[clase]);
        valor[0]=positivos-negativos;

      if (positivos==0 && menos_infinito)
        valor[0]=-2*ne;
    }
    else
      valor[0]=-2*ne;
  }
  else if (aciertos+aciertos_adicionales>granularidad){
           if (!reducir_reglas){
             valor[0]=aciertos+aciertos_adicionales-fallos-granularidad;
             if (valor[0]==0)
                valor[0]=-ne;
           }
           else {
             valor[0]=aciertos+aciertos_adicionales-fallos-granularidad;
             if (valor[0]==0)
                valor[0]=-ne;
             else if ( reglas_eliminables>=(-valor[0])+2)
                    valor[0]=aciertos+aciertos_adicionales-fallos;
                 else
                   valor[0]=-2*ne;
           }
       }
       else
           //valor[0]=aciertos-fallos;
           valor[0]=-(ne+fallos);


   valor[6]= aciertos_adicionales;
   valor[5]= aciertos-fallos;
   valor[4]= variables_estables/*+(0.00001)*unos_variables*/;
   valor[3]= (variables_irrelevantes+funciones_irrelevantes);
   //valor[2]= -negativos_global;
   valor[2]=0;
   //valor[1]= afecta;
   valor[1]= reglas_eliminables;
   /* if (reglas_eliminables>=1 and valor[0]>0){
      valor[2]=valor[0];
      valor[0]=valor[0]+((1.0*aciertos_adicionales)/ponderacion_por_clase[clase]);
      valor[4]=reglas_otra_clase;
   }*/

   /*if ( (reglas_eliminables>=1 and reglas_eliminables<2) and (valor[0]==-ne or valor[0]==0) ){
     if (!discriminante){
       valor[1]=1;
       valor[2]=balance;
       valor[3]=-ganancia_en_variables;
       //valor[4]=0;
       valor[4]=-ganancia_funciones_en_variables;
       valor[5]=diferencia_estables;
       valor[6]=max_porcen;
     }
     else {
       valor[1]=1;
       //valor[3]=0;
       if (aciertos_adicionales>0)
         valor[2]=0.1+aciertos_adicionales;
       else
         if (balance>0 and ganancia_en_variables>0){
           valor[2]=1;
         }
         else {
           valor[2]=0;
         }
       valor[3]=ganancia_en_variables;
       valor[4]=ganancia_funciones_en_variables;
       valor[5]=diferencia_estables;
       valor[6]=max_porcen;
       //valor[2]=balance;
     }
   }*/
   delete [] lista_reglas_eliminables;
  }
  else {
    valor[0]=valor[1]=valor[2]=valor[3]=valor[4]=valor[5]=valor[6]=-2*ne;
    peso=0;
  }


  if (nrule>0){
    delete [] basica_esta_regla;
    delete [] errores_esta_regla;
  }

  delete [] ejemplos_afectado_por_esta_regla;
  delete [] nr2;

  return peso;
}








double Aciertos(const ruleset &R, const VectorVar &V, example_set &E, double *agregado, double *peso_agregado, string nomfich, VectorVar &F)
{

    int n=E.N_Examples();
    int clase,kk;
    vectordouble w;
    int *marcar = new int[n];
    double *g_pos_grado = new double[n];
    double *g_neg_grado = new double[n];
    double *peso_pos = new double[n];
    double *peso_neg = new double[n];
    double bien=0, mal=0, nosesabe=0, grado;
    int conse=V.Consecuente();
    int *regla_disparada = new int[n];
    int *aislada = new int[n];
    int *bien_clas=new int[R.N_rule()];
    int *mal_clas=new int[R.N_rule()];
    int *unico_rule = new int [R.N_rule()];
    int clase_ejemplo;

    for (int i=0; i<R.N_rule(); i++)
    {
        bien_clas[i]=0;
        mal_clas[i]=0;
        unico_rule[i]=0;
    }


    for (int i=0; i<n; i++)
    {
        w=E.Data(i);
        clase_ejemplo=(int) (w.At(conse));
        clase=R.InferenceC(w,clase_ejemplo,g_pos_grado[i],peso_pos[i],g_neg_grado[i],peso_neg[i],regla_disparada[i],aislada[i]);


        if (clase==-1)
        {
            nosesabe++;
            //agregado[i]=0;
            //peso_agregado[i]=0;
            marcar[i]=-1;
        }
        else if (clase==w.At(conse))
        {
            //else if (V.Adaptacion(w.At(conse),conse,clase)>0){
            bien++;
            marcar[i]=1;
            //agregado[i]=grado;
            //peso_agregado[i]=R.Get_Weight(regla_disparada[i]);
            bien_clas[regla_disparada[i]]++;
            if (aislada[i]==1)
              unico_rule[regla_disparada[i]]++;
        }
        else
        {
            mal++;
            marcar[i]=0;
            //agregado[i]=-grado;
            //peso_agregado[i]=R.Get_Weight(regla_disparada[i]);
            mal_clas[regla_disparada[i]]++;
        }
    }

    E.Mark(marcar,n, g_pos_grado, peso_pos, g_neg_grado, peso_neg, regla_disparada, aislada);
    //E.UnMark(desmarcar,tama2,gdesmarcar);
    ofstream out(nomfich.c_str(),ios::app);

    out <<"Aciertos: " << bien << endl;
    out <<"Errores:  " << mal << endl;
    out <<"No clasificados: " << nosesabe << endl;

    int* Necesaria_Rule = new int [R.N_rule()];
    for (int i=0; i<R.N_rule(); i++)
      Necesaria_Rule[i]=0;

    for (int i=0; i<n; i++){
      if (regla_disparada[i]!=-1 and aislada[i]==1)
        Necesaria_Rule[regla_disparada[i]]++;
    }


    for (int i=0; i<R.N_rule(); i++)
    {
        out << "\tRule "<<i <<"("<<R.SPrint_Consequent(i) <<"): " << bien_clas[i] << " / " << mal_clas[i] << " (basico: " << unico_rule[i] <<")" << endl;
        cout << "\tRule "<<i <<"("<<R.SPrint_Consequent(i) <<"): " << bien_clas[i] << " / " << mal_clas[i] << " (basico: " << unico_rule[i] <<")" << endl;
        //cout << "\tRule "<<i <<": " <<Necesaria_Rule[i] << endl;
    }

    cout << "Numero de reglas: " << R.N_rule() << endl;

    out.close();

    //cin.get();

    delete [] Necesaria_Rule;

    delete []marcar;
    delete []g_pos_grado;
    delete []g_neg_grado;
    delete []peso_pos;
    delete []peso_neg;
    delete []bien_clas;
    delete []regla_disparada;
    delete []aislada;
    delete []unico_rule;
    delete []mal_clas;

    return bien/n;

}







example_set RetenerEjemplos(ruleset &R, const VectorVar &V, example_set &E, VectorVar &F)
{
    int n=E.N_Examples();
    bool *tomar = new bool[R.N_rule()];
    vectordouble w;
    vector<int> s;

    for (int i=0; i<R.N_rule(); i++)
      tomar[i]=true;



    for (int i=0; i<n; i++)
    {
      if (E.Rule_Used(i)==-1)
	s.push_back(i);
      else if (tomar[E.Rule_Used(i)]){
	s.push_back(i);
	tomar[E.Rule_Used(i)]=false;
      }
      else {
	tomar[E.Rule_Used(i)]=true;
      }
    }


    example_set E_result(s.size(),V.N_Antecedente()+1);
    for (int i=0; i<s.size(); i++){
      w = E.Data(s[i]);
      E_result.AddExample(w,i);
    }

    return E_result;

}

double AciertosLigero(ruleset &R, const VectorVar &V, example_set &E, VectorVar &F)
{
    int n=E.N_Examples();
    cout << "@" << n << endl << endl;
    int clase,kk;
    vectordouble w;

    double g_pos_grado;
    double g_neg_grado;
    double peso_pos;
    double peso_neg;
    int tama2=0;
    double bien=0, mal=0, nosesabe=0, grado;
    int conse=V.Consecuente();
    int regla_disparada;
    int aislada;
    int *bien_clas=new int[R.N_rule()];
    int *mal_clas=new int[R.N_rule()];
    int *unico_rule = new int [R.N_rule()];
    int clase_ejemplo;

    for (int i=0; i<R.N_rule(); i++)
    {
        bien_clas[i]=0;
        mal_clas[i]=0;
        unico_rule[i]=0;
    }



    for (int i=0; i<n; i++)
    {
        w=E.Data(i);
        clase_ejemplo=(int) (w.At(conse));
				deleteLine();
				cout << i << "\t" << clase_ejemplo << "\t" << endl;
				clase=R.InferenceC(w,clase_ejemplo,g_pos_grado,peso_pos,g_neg_grado,peso_neg,regla_disparada,aislada);
        if (clase==-1)
        {
            nosesabe++;
        }
        else if (clase==w.At(conse))
        {
            //else if (V.Adaptacion(w.At(conse),conse,clase)>0){
            bien++;
            bien_clas[regla_disparada]++;
            if (aislada==1)
              unico_rule[regla_disparada]++;
        }
        else
        {
            mal++;
            mal_clas[regla_disparada]++;
        }
				deleteLine();
				cout << i << "\t" << clase_ejemplo << "\t" << clase << endl;
    }



    delete []bien_clas;
    delete []mal_clas;
    delete []unico_rule;

    return bien/n;

}



double Aciertos(ruleset &R, const VectorVar &V, example_set &E, VectorVar &F)
{
    int n=E.N_Examples();
    cout << "@" << n << endl << endl;
    int clase,kk;
    vectordouble w;
    int *marcar = new int[n], tama1=0;
    double *g_pos_grado = new double[n];
    double *g_neg_grado = new double[n];
    double *peso_pos = new double[n];
    double *peso_neg = new double[n];
    int tama2=0;
    double bien=0, mal=0, nosesabe=0, grado;
    int conse=V.Consecuente();
    int *regla_disparada = new int[n];
    int *aislada = new int[n];
    int *bien_clas=new int[R.N_rule()];
    int *mal_clas=new int[R.N_rule()];
    int *unico_rule = new int [R.N_rule()];
    int clase_ejemplo;

    for (int i=0; i<R.N_rule(); i++)
    {
        bien_clas[i]=0;
        mal_clas[i]=0;
        unico_rule[i]=0;
    }



    for (int i=0; i<n; i++)
    {
        w=E.Data(i);
        clase_ejemplo=(int) (w.At(conse));
        clase=R.InferenceC(w,clase_ejemplo,g_pos_grado[i],peso_pos[i],g_neg_grado[i],peso_neg[i],regla_disparada[i],aislada[i]);
        if (clase==-1)
        {
            nosesabe++;
            marcar[i]=-1;
        }
        else if (clase==w.At(conse))
        {
            //else if (V.Adaptacion(w.At(conse),conse,clase)>0){
            bien++;
            marcar[i]=1;
            tama1++;
            bien_clas[regla_disparada[i]]++;
            if (aislada[i]==1)
              unico_rule[regla_disparada[i]]++;
        }
        else
        {
            mal++;
            marcar[i]=0;
            tama2++;
            mal_clas[regla_disparada[i]]++;
        }
				deleteLine();
				cout << i << "\t" << clase_ejemplo << "\t" << clase << endl;
    }

    E.Mark(marcar,n, g_pos_grado, peso_pos, g_neg_grado, peso_neg, regla_disparada, aislada);
    //E.UnMark(desmarcar,tama2,gdesmarcar);

    /*cout <<"Aciertos: " << bien << endl;
    cout <<"Errores:  " << mal << endl;
    cout <<"No clasificados: " << nosesabe << endl;

    for (int i=0; i<R.N_rule(); i++)
    {
        cout << "\tRule "<<i <<"("<<R.SPrint_Consequent(i) <<"): " << bien_clas[i] << " / " << mal_clas[i] << " (basico: " << unico_rule[i] <<")" << endl;
    }*/



    delete []marcar;
    delete []g_pos_grado;
    delete []g_neg_grado;
    delete []peso_pos;
    delete []peso_neg;
    delete []bien_clas;
    delete []mal_clas;
    delete []unico_rule;
    delete []regla_disparada;
    delete []aislada;

    return bien/n;

}


double Aciertos(ruleset &R, const VectorVar &V, example_set &E, string nomfich, VectorVar &F)
{
    int n=E.N_Examples();
    int clase,kk;
    vectordouble w;
    int *marcar = new int[n];
    double *g_pos_grado = new double[n];
    double *g_neg_grado = new double[n];
    double *peso_pos = new double[n];
    double *peso_neg = new double[n];
    double bien=0, mal=0, nosesabe=0, grado;
    int conse=V.Consecuente();
    int *regla_disparada = new int[n];
    int *aislada = new int[n];
    int *bien_clas=new int[R.N_rule()];
    int *mal_clas=new int[R.N_rule()];
    int clase_ejemplo;

    for (int i=0; i<R.N_rule(); i++)
    {
        bien_clas[i]=0;
        mal_clas[i]=0;
    }



    for (int i=0; i<n; i++)
    {
        w=E.Data(i);
        clase_ejemplo=(int) (w.At(conse));
        clase=R.InferenceC(w,clase_ejemplo,g_pos_grado[i],peso_pos[i],g_neg_grado[i],peso_neg[i],regla_disparada[i],aislada[i]);
        if (clase==-1)
        {
            nosesabe++;
            marcar[i]=-1;
        }
        else if (clase==w.At(conse))
        {
            //else if (V.Adaptacion(w.At(conse),conse,clase)>0){
            bien++;
            marcar[i]=1;
             bien_clas[regla_disparada[i]]++;
        }
        else
        {
            mal++;
            marcar[i]=0;
            mal_clas[regla_disparada[i]]++;
        }
    }

    E.Mark(marcar,n, g_pos_grado, peso_pos, g_neg_grado, peso_neg, regla_disparada, aislada);
    //E.UnMark(desmarcar,tama2,gdesmarcar);
    ofstream out(nomfich.c_str(),ios::app);

    out <<"Aciertos: " << bien << endl;
    out <<"Errores:  " << mal << endl;
    out <<"No clasificados: " << nosesabe << endl;

    for (int i=0; i<R.N_rule(); i++)
    {
        out << "\tRule "<<i <<"("<<R.SPrint_Consequent(i) <<"): " << bien_clas[i] << " / " << mal_clas[i] << endl;
    }

    out.close();


    delete []marcar;
    delete []g_pos_grado;
    delete []g_neg_grado;
    delete []peso_pos;
    delete []peso_neg;
    delete []bien_clas;
    delete []mal_clas;
    delete []regla_disparada;
    delete []aislada;


    return bien/n;

}









double Aciertos_Lista_Reglas_Eliminables(ruleset &R, const VectorVar &V, example_set &E, VectorVar &F, int*lista, int &tama, int *l_aciertos, int &tama_l)
{
    int n=E.N_Examples();
    int clase,kk;
    vectordouble w;
    int *marcar = new int[n], tama1=0;
    double *g_pos_grado = new double[n];
    double *g_neg_grado = new double[n];
    double *peso_pos = new double[n];
    double *peso_neg = new double[n];
    int tama2=0;
    double bien=0, mal=0, nosesabe=0, grado;
    int conse=V.Consecuente();
    int *regla_disparada = new int[n];
    int *aislada = new int[n];
    int *bien_clas=new int[R.N_rule()];
    int *mal_clas=new int[R.N_rule()];
    int *unico_rule = new int [R.N_rule()];
    bool *usada = new bool[R.N_rule()];
    int clase_ejemplo;

    for (int i=0; i<R.N_rule(); i++)
    {
        bien_clas[i]=0;
        mal_clas[i]=0;
        unico_rule[i]=0;
	usada[i]=false;
    }


    tama_l=0;
    for (int i=0; i<n; i++)
    {
	l_aciertos[i]=0;
        w=E.Data(i);
        clase_ejemplo=(int) (w.At(conse));
        clase=R.InferenceC(w,clase_ejemplo,g_pos_grado[i],peso_pos[i],g_neg_grado[i],peso_neg[i],regla_disparada[i],aislada[i]);
        if (clase==-1)
        {
            nosesabe++;
            marcar[i]=-1;
        }
        else if (clase==w.At(conse))
        {
            //else if (V.Adaptacion(w.At(conse),conse,clase)>0){
            bien++;
            marcar[i]=1;
            tama1++;
            bien_clas[regla_disparada[i]]++;
	    l_aciertos[tama_l]=i;
	    tama_l++;
	    usada[regla_disparada[i]]=true;
            if (aislada[i]==1)
              unico_rule[regla_disparada[i]]++;
        }
        else
        {
            mal++;
            marcar[i]=0;
            tama2++;
            mal_clas[regla_disparada[i]]++;
	    usada[regla_disparada[i]]=true;
        }
    }

    /*E.Mark(marcar,n, g_pos_grado, peso_pos, g_neg_grado, peso_neg, regla_disparada, aislada);
    E.UnMark(desmarcar,tama2,gdesmarcar);

    cout <<"Aciertos: " << bien << endl;
    cout <<"Errores:  " << mal << endl;
    cout <<"No clasificados: " << nosesabe << endl;

    for (int i=0; i<R.N_rule(); i++)
    {
        cout << "\tRule "<<i <<"("<<R.SPrint_Consequent(i) <<"): " << bien_clas[i] << " / " << mal_clas[i] << " (basico: " << unico_rule[i] <<")" << endl;
    }*/

    tama=0;
    for (int i=0; i<R.N_rule(); i++){
      if (usada[i]){
        lista[tama]=i;
        tama++;
      }
    }


    delete []marcar;
    delete []g_pos_grado;
    delete []g_neg_grado;
    delete []peso_pos;
    delete []peso_neg;
    delete []bien_clas;
    delete []mal_clas;
    delete []unico_rule;
    delete []regla_disparada;
    delete []aislada;
    delete []usada;

    return bien/n;

}















double Inferencia(ruleset &R, const VectorVar &V, example_set &E, VectorVar &F)
{
    int n=E.N_Examples();
    int clase, regla_disparada;
    vectordouble w;
    int conse=V.Consecuente();
    double nosesabe=0, bien=0, mal=0, grado;
    for (int i=0; i<n; i++)
    {
        w=E.Data(i);
        clase=R.InferenceC(w,grado,regla_disparada);
        if (clase==-1)
        {
            nosesabe++;
        }
        else if (clase==w.At(conse))
        {
            bien++;
        }
        else
        {
            mal++;
        }
    }
    return bien/n;
}





double Inferencia2BasesReglas(ruleset &R, ruleset &R_Descrip, const VectorVar &V, example_set &E, int *fallosR, int* fallosR_Descrip, int *disparosR, int *disparosR_Descrip)
{
    int n=E.N_Examples();
    int clase, clase1, regla_disparada1, clase2, regla_disparada2, discrepancias=0, discrepancias_reglas=0;

    if (R.N_rule()!=R_Descrip.N_rule())
      cout << "\nError: No coinciden las bases de reglas en el numero de reglas. Descriptivas: " << R_Descrip.N_rule() << " Discriminantes: " << R.N_rule() << endl;

    for (int i=0; i<R.N_rule(); i++){
      fallosR[i]=0;
      fallosR_Descrip[i]=0;
      disparosR[i]=0;
      disparosR_Descrip[i]=0;
      //Pinta2Reglas(R_Descrip,R,V,i);

    }
    vectordouble w;
    int conse=V.Consecuente();
    double nosesabe=0, bien=0, mal=0, grado1, grado2, fallos1=0, fallos2=0, n1=0, n2=0;
    for (int i=0; i<n; i++)
    {
        w=E.Data(i);
        clase1=R.InferenceC(w,grado1,regla_disparada1);
        clase2=R_Descrip.InferenceC(w,grado2,regla_disparada2);

	if (regla_disparada1!=-1)
	  disparosR[regla_disparada1]++;

	if (regla_disparada2!=-1)
	  disparosR_Descrip[regla_disparada2]++;

	if (regla_disparada1!=regla_disparada2)
	  discrepancias_reglas++;

	if (clase1!=-1 and clase1!=w.At(conse))
	  fallosR[regla_disparada1]++;

	if (clase2!=-1 and clase2!=w.At(conse))
	  fallosR_Descrip[regla_disparada2]++;

	if (clase1!=-1)
	  n1++;

	if (clase2!=-1)
	  n2++;

	if (clase1!=-1 and clase1!=w.At(conse))
	  fallos1++;

	if (clase2!=-1 and clase2!=w.At(conse))
	  fallos2++;


        if (clase1==clase2){
            clase=clase1;
	}
        else
        {
            discrepancias++;
            if (clase2==-1)
                clase = clase1;
            else if (clase1 == -1)
                clase = clase2;
            else
            {
                if (grado1>grado2)
                    clase = clase1;
                else
                    clase = clase2;

            }
        }

        if (clase==-1)
        {
            nosesabe++;
        }
        else if (clase==w.At(conse))
        {
            bien++;
        }
        else
        {
	    mal++;
        }
    }
    /*cout << "Discrepancias clases= " << discrepancias << endl;
    cout << "Discrepancias reglas= " << discrepancias_reglas << endl;
    if (n1>0)
      cout << "Fallos R1= " << (1.0*fallos1)/n1 << " ("<< fallos1 << "/" << n1 <<")\n";

    if (n2>0)
      cout << "Fallos R2= " << (1.0*fallos2)/n2 << " ("<< fallos2 << "/" << n2 <<")\n";

    cout << "Acierto = " << bien/n << " ("<< bien << "/" << n <<")\n";


    for (int i=0; i<R.N_rule(); i++){
      cout << "Fallos regla " << i << ": " << fallosR[i] << " | " << fallosR_Descrip[i] << " || DisR: " << disparosR[i] << " | " << disparosR_Descrip[i] << endl;
    }*/

    //char ch;
    //cin >> ch;


    return bien/n;
}



double Inferencia2BasesReglas_2aprox(ruleset &R, ruleset &R_Descrip, const VectorVar &V, example_set &E, int *fallosR, int* fallosR_Descrip, int *disparosR, int *disparosR_Descrip)
{
    int n=E.N_Examples();
    int clase, clase1, regla_disparada1, clase2, regla_disparada2, discrepancias=0, discrepancias_reglas=0;

    if (R.N_rule()!=R_Descrip.N_rule())
      cout << "\nError: No coinciden las bases de reglas en el numero de reglas. Descriptivas: " << R_Descrip.N_rule() << " Discriminantes: " << R.N_rule() << endl;

    for (int i=0; i<R.N_rule(); i++){
      fallosR[i]=0;
      fallosR_Descrip[i]=0;
      disparosR[i]=0;
      disparosR_Descrip[i]=0;
      //Pinta2Reglas(R_Descrip,R,V,i);

    }
    vectordouble w;
    int conse=V.Consecuente();
    double nosesabe=0, bien=0, mal=0, grado1, grado2, fallos1=0, fallos2=0, n1=0, n2=0;
    for (int i=0; i<n; i++)
    {
        w=E.Data(i);
        clase1=R.InferenceC(w,grado1,regla_disparada1);
        clase2=R_Descrip.InferenceC(w,grado2,regla_disparada2);

	if (regla_disparada1!=-1)
	  disparosR[regla_disparada1]++;

	if (regla_disparada2!=-1)
	  disparosR_Descrip[regla_disparada2]++;

	if (regla_disparada1!=regla_disparada2)
	  discrepancias_reglas++;

	if (clase1!=-1 and clase1!=w.At(conse))
	  fallosR[regla_disparada1]++;

	if (clase2!=-1 and clase2!=w.At(conse))
	  fallosR_Descrip[regla_disparada2]++;

	if (clase1!=-1)
	  n1++;

	if (clase2!=-1)
	  n2++;

	if (clase1!=-1 and clase1!=w.At(conse))
	  fallos1++;

	if (clase2!=-1 and clase2!=w.At(conse))
	  fallos2++;


        if (clase1==clase2){
            clase=clase1;
	}
        else
        {
            discrepancias++;
            if (clase2==-1)
                clase = clase1;
            else if (clase1 == -1)
                clase = clase2;
            else
            {
                if (grado1>grado2)
                    clase = clase1;
                else
                    clase = clase2;

            }
        }

        if (clase==-1)
        {
            nosesabe++;
        }
        else if (clase==w.At(conse))
        {
            bien++;
        }
        else
        {
	    mal++;
        }
    }
    /*cout << "Discrepancias clases= " << discrepancias << endl;
    cout << "Discrepancias reglas= " << discrepancias_reglas << endl;
    if (n1>0)
      cout << "Fallos R1= " << (1.0*fallos1)/n1 << " ("<< fallos1 << "/" << n1 <<")\n";

    if (n2>0)
      cout << "Fallos R2= " << (1.0*fallos2)/n2 << " ("<< fallos2 << "/" << n2 <<")\n";

    cout << "Acierto = " << bien/n << " ("<< bien << "/" << n <<")\n";


    for (int i=0; i<R.N_rule(); i++){
      cout << "Fallos regla " << i << ": " << fallosR[i] << " | " << fallosR_Descrip[i] << " || DisR: " << disparosR[i] << " | " << disparosR_Descrip[i] << endl;
    }*/

    //char ch;
    //cin >> ch;


    return bien/n;
}





double Filtrar_Reglas(ruleset &R, const VectorVar &V, example_set &E, VectorVar &F)
{
    double mejor_acierto, acierto, peso;
    string regla_escrita;

    mejor_acierto = Inferencia(R,V,E,F);
    if (R.N_rule()>2)
    {
        cout << "La Regla:\n";
        cout << R.SPrint(R.N_rule()-1)<< endl;
        cout << "################## POR ###############\n";

        int j=0;
        genetcode *aux;
        //cout <<"\t\tMejor Aciertos: "<< mejor_acierto << endl;
        while (j<R.N_rule())
        {
            regla_escrita = R.SPrint(j);
            aux = R.Extract(j,peso);
            acierto = Inferencia(R,V,E,F);
            //cout << "\t\tSi quito la regla " << j <<": " <<acierto << endl;
            if (acierto>=mejor_acierto)
            {
                delete aux;
                mejor_acierto=acierto;
                //cout << "Eliminada regla " << j << " acierto= " << mejor_acierto << endl;
                //cout << "Numero de reglas: " << R.N_rule() << endl;
                j=0;
                cout << regla_escrita << endl;
            }
            else
            {
                R.Insert(j,aux,peso);
                //acierto = Inferencia(R,V,E);
                //cout << "Numero de reglas: " << R.N_rule() << endl;
                j++;
            }
        }
        //cin.get();
        //sleep(2);
    }

    return mejor_acierto;
}


double Filtrar_Reglas_New(ruleset &R, ruleset &R_dis, ruleset &R_mix, int *reglas_a_eliminar,
			  int n_reglas_a_eliminar, bool Specific_Rule)
{
    double mejor_acierto, acierto, peso;
    string regla_escrita;

    /*if (!Specific_Rule){*/

      for (int i=n_reglas_a_eliminar-1; i>=0; i--){
	genetcode *aux;
	//regla_escrita = R.SPrint(reglas_a_eliminar[i]);
	//cout << regla_escrita << endl;
	aux = R.Extract(reglas_a_eliminar[i],peso);
	delete aux;
	aux = R_dis.Extract(reglas_a_eliminar[i],peso);
	delete aux;
	aux = R_mix.Extract(reglas_a_eliminar[i],peso);
	delete aux;
      }
    /*}
    else {
      genetcode *aux;
      aux = R.Extract(reglas_a_eliminar[0],peso);
      delete aux;
      R_dis.Remove();
      aux = R_dis.Extract(reglas_a_eliminar[0],peso);
      R_dis.Add(aux[0],peso);
    }*/


    return 0;
}




double Filtrar_Reglas(ruleset &R, const VectorVar &V, example_set &E, ruleset &R_Descrip, VectorVar &F)
{
    double mejor_acierto, acierto, peso, peso2;

    mejor_acierto = Inferencia(R,V,E,F);
    if (R.N_rule()>2)
    {
        int j=0;
        genetcode *aux, *aux2;
        while (j<R.N_rule())
        {
            aux = R.Extract(j,peso);
            acierto = Inferencia(R,V,E,F);
            if (acierto>=mejor_acierto)
            {
                delete aux;
                mejor_acierto=acierto;
                aux2=R_Descrip.Extract(j,peso2);
                delete aux2;
                cout << "Eliminada regla " << j << " acierto= " << mejor_acierto << endl;
                cout << "Numero de reglas: " << R.N_rule() << endl;
                j=0;
            }
            else
            {
                R.Insert(j,aux,peso);
                //acierto = Inferencia(R,V,E);
                //cout << "Numero de reglas: " << R.N_rule() << endl;
                j++;
            }
        }
    }

    return mejor_acierto;
}









void PintaIndividuo(genetcode GC, const VectorVar &V)
{

}



// devuelve una matriz de nv filas (n. variables antecedente) y
// m columnas (numero de clases del consecuente +1) donde en la
// posicion 0 se encuentra la medida general sobre todas las clases.
double **ReservaParaMedidasDeInformacion(const VectorVar &V)
{
cout << "Med " << V.N_Antecedente() << " " << V.SizeDomain(V.Consecuente()) + 1 << endl;
    int n=V.TotalVariables(), conse=V.Consecuente();
    int tam_dom_conse=V.SizeDomain(conse);
    int nv, m=tam_dom_conse+1;


    nv =V.N_Antecedente();
    // Reservo una matriz de tamano nvxm
    double **matriz= new double *[nv];
    for (int i=0; i<nv; i++) {
        matriz[i]=new double [m];
//	for (int j = 0; j < m; j++)
//		matriz[i][j] = 0;
    }

    return matriz;
}

// Libera la matriz de medidas de informacion
void LiberaMedidasDeInformacion(double **&matriz, int n)
{
    for (int i=0; i<n; i++)
        delete [] matriz[i];

    delete [] matriz;
}



void CalcularProbabilidadSimpleVariable(const VectorVar &V, const example_set &E, int variable, int n_casos, double *p)
{
    int ne=0;
    double suma=0;
    double *aux=new double[n_casos];

    for (int j=0; j<n_casos; j++)
        p[j]=0;

    for (int i=0; i<E.N_Examples(); i++)
    {
        if (!E.Is_Covered(i))
        {
            suma=0;
            for (int j=0; j<n_casos; j++)
            {
                aux[j]=V.Domain(variable).Adaptacion(E.Data(i,variable),j);
                suma+=aux[j];
            }
            // A�dir los valores a la sumatoria
            if (suma>0)
            {
                for (int j=0; j<n_casos; j++)
                    p[j]=p[j]+(aux[j]/suma);

                ne++;
            }
        }
    }

    delete [] aux;

    if (ne>0)
    {
      	      cout << "NE: " << ne << endl;

        for (int j=0; j<n_casos; j++)
            p[j]=p[j]/ne;
    }

}


double Operacion (double x1, double x2, int op)
{
    double a;
    if (x1==-999999999 || x2==-999999999)
        return -999999999;
    switch (op)
    {
        //Operacion ==
    case 19:
        if (x1 == x2)
            return 1;
        else
            return 0;
        break;
        //Operacion OR
    case 20:
        if (x1 == 1 || x2 == 1)
            return 1;
        else
            return 0;
        break;
        //Operacion XOR
    case 21:
        if ( (x1 == 1 && x2 == 0 ) || (x1 == 0 && x2 == 1 ) )
            return 1;
        else
            return 0;
        break;
    case 22:
        if ( (x1 == 1 && x2 == 1 ))
            return 1;
        else
            return 0;
        break;

    case 1:
        if (x1>=x2)
            return 1;
        else
            return 0;
        break;

    case 2:
        return x1+x2;
        break;

        /*   case 5: return x1-x2;
                   break;
        */
    case 5: // Division
        if (x2!=0)
           return x1/x2;
        else
           return -999999999;
        break;

    case 4:
        return (x1-x2);
        break;

    case 3:
        return (x1*x2);
        break;
    case 7: //Delta+
        /*a = sin(x1);
        if (a<-1 || a >1)
            cout << "sin(" << x1 << ")=" << a << endl;*/
        return x1;
        break;
    case 8: //Delta-
        return x1;
        break;
    case 9:
        return log(fabs(x1)+1);
        break;
    case 10:
        return sin(x1)+cos(x1);
        break;
    case 11:
        return sqrt(fabs(x1));
        break;
    case 12:
        return 1.0/(fabs(x1)+1);
        break;
    case 13:
        break;
    case 15:
        return x1;
        break;
    }
}

void CalcularProbabilidadSimpleVariable_op(const VectorVar &V, const example_set &E, int v1, int v2,
        domain_t &atrib, int op, double *p)
{
    int ne=0;
    int n_casos = atrib.Size();
    double suma=0;
    double *aux=new double[n_casos];
    double example;
    bool error;

    for (int j=0; j<n_casos; j++)
        p[j]=0;

    double min, max;
    if (op==6 || op==13)
    {
        Max_Min(E,V,v1,v2,-1,min,max);
    }

    for (int i=0; i<E.N_Examples(); i++)
    {
        if (!E.Is_Covered(i))
        {
            suma=0;
            error=false;

            if (op!=6 && op!=13){
              if (op==5 and E.Data(i,v2)==0)
                error=true;
              else
                example = Operacion(E.Data(i,v1),E.Data(i,v2), op);
            }
            else if (op==6)
                example = (E.Data(i,v1)-min) / (max-min);
            else
                example = 2* (E.Data(i,v1)-min) / (max-min)-1;

            for (int j=0; j<n_casos and !error; j++)
            {
                aux[j]=atrib.Adaptacion(example,j);
                suma+=aux[j];
            }
            // A�dir los valores a la sumatoria
            if (suma>0)
            {
                for (int j=0; j<n_casos; j++)
                    p[j]=p[j]+(aux[j]/suma);

                ne++;
            }
        }
    }

    delete [] aux;

    if (ne>0)
    {
        for (int j=0; j<n_casos; j++)
            p[j]=p[j]/ne;
    }


    //for (int j=0; j<n_casos; j++)
    //cout << "p(" << j << ")= " << p[j] << endl;

}



void CalcularProbabilidadConjunta2Variables(const VectorVar &V, const example_set &E, int var1, int var2,
        int n_casos1, int n_casos2, double **m)
{

    int ne=0;
    double suma;

    double **aux = new double*[n_casos1];
    for (int j=0; j<n_casos1; j++)
        aux[j] = new double[n_casos2];

    for (int j=0; j<n_casos1; j++)
        for (int k=0; k<n_casos2; k++)
            m[j][k]=0;


    for (int i=0; i<E.N_Examples(); i++)
    {
        if (!E.Is_Covered(i))
        {
            suma=0;
            for (int j=0; j<n_casos1; j++)
                for (int k=0; k<n_casos2; k++)
                {
                    aux[j][k]=V.Domain(var1).Adaptacion(E.Data(i,var1),j)*V.Domain(var2).Adaptacion(E.Data(i,var2),k);
                    suma=suma+aux[j][k];
                }

            if (suma>0)
            {
                for (int j=0; j<n_casos1; j++)
                    for (int k=0; k<n_casos2; k++)
                        m[j][k]=m[j][k]+(aux[j][k]/suma);

            }

           ne++;
        }

    }

    if (ne>0)
        for (int j=0; j<n_casos1; j++)
            for (int k=0; k<n_casos2; k++)
                m[j][k]=m[j][k]/ne;

    /*for (int j=0; j<n_casos1; j++)
      for (int k=0; k<n_casos2; k++)
        cout << "p(x=" << j <<", y=" << k <<")= " << m[j][k] << endl; */


    for (int j=0; j<n_casos1; j++)
        delete [] aux[j];

    delete [] aux;
}




void CalcularProbabilidadConjunta3Variables(const VectorVar &V, const example_set &E, int var1, int var2, int var3,
        int n_casos1, int n_casos2, int n_casos3, double ***m)
{

    int ne=0;
    double suma;

    double ***aux = new double**[n_casos1];
    for (int j=0; j<n_casos1; j++)
    {
        aux[j] = new double*[n_casos2];
        for (int k=0; k<n_casos2; k++)
            aux[j][k]= new double[n_casos3];
    }


    for (int j=0; j<n_casos1; j++)
        for (int k=0; k<n_casos2; k++)
            for (int t=0; t<n_casos3; t++)
                m[j][k][t]=0;


    for (int i=0; i<E.N_Examples(); i++)
    {
        if (!E.Is_Covered(i))
        {
            suma=0;
            for (int j=0; j<n_casos1; j++)
                for (int k=0; k<n_casos2; k++)
                    for (int t=0; t<n_casos3; t++)
                    {
                        aux[j][k][t]=V.Domain(var1).Adaptacion(E.Data(i,var1),j)*V.Domain(var2).Adaptacion(E.Data(i,var2),k)*V.Domain(var3).Adaptacion(E.Data(i,var3),t);
                        suma=suma+aux[j][k][t];
                    }

            if (suma>0)
            {
                for (int j=0; j<n_casos1; j++)
                    for (int k=0; k<n_casos2; k++)
                        for (int t=0; t<n_casos3; t++)
                            m[j][k][t]=m[j][k][t]+(aux[j][k][t]/suma);

                ne++;
            }
        }

    }

    if (ne>0)
        for (int j=0; j<n_casos1; j++)
            for (int k=0; k<n_casos2; k++)
                for (int t=0; t<n_casos3; t++)
                    m[j][k][t]=m[j][k][t]/ne;

    /*for (int j=0; j<n_casos1; j++)
      for (int k=0; k<n_casos2; k++)
        cout << "p(x=" << j <<", y=" << k <<")= " << m[j][k] << endl; */


    for (int j=0; j<n_casos1; j++)
    {
        for (int k=0; k<n_casos2; k++)
            delete [] aux[j][k];
        delete [] aux[j];
    }
    delete [] aux;
}




double CalcularInformacionCondicionada2variablesRespectoClase (const VectorVar &V, const example_set &E, int casos1, int casos2,
        int casosclase, double ***Conjunta123, double **Conjunta13,
        double **Conjunta23, double *simple3, double *I)
{
    double Inf_global=0, Inf_parcial;

    for (int z=0; z<casosclase; z++)
    {
        Inf_parcial=0;
        for (int y=0; y<casos2; y++)
        {
            for (int x=0; x<casos1; x++)
            {
                if (Conjunta123[x][y][z]!=0 && Conjunta13[x][z]*Conjunta23[y][z]!=0)
                    Inf_parcial+=Conjunta123[x][y][z]*log( (simple3[z]*Conjunta123[x][y][z]) / (Conjunta13[x][z]*Conjunta23[y][z]));
            }
        }
        I[z]=Inf_parcial;
        Inf_global+=Inf_parcial;
    }
    return Inf_global;
}





void CalcularProbabilidadConjunta2Variables_op(const VectorVar &V, const example_set &E, int var1, int var2,
        domain_t &atrib, int vy, int n_casos2, int op, double **m)
{

    int ne=0;
    int n_casos1=atrib.Size();
    double example;
    double suma;
    bool error=false;

    double **aux = new double*[n_casos1];
    for (int j=0; j<n_casos1; j++)
        aux[j] = new double[n_casos2];

    for (int j=0; j<n_casos1; j++)
        for (int k=0; k<n_casos2; k++)
            m[j][k]=0;


    double min, max;
    if (op==6 || op==13)
    {
        Max_Min(E,V,var1,var2,-1,min,max);
    }


    for (int i=0; i<E.N_Examples(); i++)
    {
        if (!E.Is_Covered(i))
        {
            suma=0;
            error=false;
            if (op!=6  &&  op!=13){
                double dato1, dato2;
                dato1=E.Data(i,var1);
                dato2=E.Data(i,var2);
                if (op==5 and dato2==0)
                  error=true;
                else
                  example = Operacion(dato1,dato2, op);

            }
            else if (op==6)
                example = (E.Data(i,var1)-min) / (max-min);
            else
                example = 2* ((E.Data(i,var1)-min) / (max-min))-1;

            for (int j=0; j<n_casos1 and !error; j++)
                for (int k=0; k<n_casos2; k++)
                {
                    aux[j][k]=atrib.Adaptacion(example,j)*V.Domain(vy).Adaptacion(E.Data(i,vy),k);
                    suma=suma+aux[j][k];
                }

            if (suma>0)
            {
                for (int j=0; j<n_casos1; j++)
                    for (int k=0; k<n_casos2; k++)
                        m[j][k]=m[j][k]+(aux[j][k]/suma);

            }
          ne++;

        }

    }


    if (ne>0)
        for (int j=0; j<n_casos1; j++)
            for (int k=0; k<n_casos2; k++)
                m[j][k]=m[j][k]/ne;


    for (int j=0; j<n_casos1; j++)
        delete [] aux[j];

    delete [] aux;
}


void CalcularProbabilidadConjunta2Variables_op_pinta(const VectorVar &V, const example_set &E, int var1, int var2,
        domain_t &atrib, int vy, int n_casos2, int op, double **m)
{

    int ne=0;
    int n_casos1=atrib.Size();
    double example;
    double suma;
    bool error=false;

    double **aux = new double*[n_casos1];
    for (int j=0; j<n_casos1; j++)
        aux[j] = new double[n_casos2];

    for (int j=0; j<n_casos1; j++)
        for (int k=0; k<n_casos2; k++)
            m[j][k]=0;


    double min, max;
    if (op==6 || op==13)
    {
        Max_Min(E,V,var1,var2,-1,min,max);
    }


    for (int i=0; i<E.N_Examples(); i++)
    {
        if (!E.Is_Covered(i))
        {
            suma=0;
            error=false;
            if (op!=6  &&  op!=13){
                double dato1, dato2;
                dato1=E.Data(i,var1);
                dato2=E.Data(i,var2);
                if (op==5 and dato2==0)
                  error=true;
                else
                  example = Operacion(dato1,dato2, op);

            }
            else if (op==6)
                example = (E.Data(i,var1)-min) / (max-min);
            else
                example = 2* ((E.Data(i,var1)-min) / (max-min))-1;

            for (int j=0; j<n_casos1 and !error; j++)
                for (int k=0; k<n_casos2; k++)
                {
                    aux[j][k]=atrib.Adaptacion(example,j)*V.Domain(vy).Adaptacion(E.Data(i,vy),k);
                    suma=suma+aux[j][k];
                }

            if (suma>0)
            {
                for (int j=0; j<n_casos1; j++)
                    for (int k=0; k<n_casos2; k++)
                        m[j][k]=m[j][k]+(aux[j][k]/suma);

            }
          ne++;
        }

    }


    if (ne>0)
        for (int j=0; j<n_casos1; j++)
            for (int k=0; k<n_casos2; k++)
                m[j][k]=m[j][k]/ne;




    cout << endl;
    for (int j=0; j<n_casos1; j++){
      cout << j;
      for (int k=0; k<n_casos2; k++){
        cout <<"\t"<< m[j][k] <<"| ";
      }
      cout << endl;
    }
    cout << endl;



    for (int j=0; j<n_casos1; j++)
        delete [] aux[j];

    delete [] aux;
}


// Calcula las medidas de informacion
void MedidasDeInformacion(const VectorVar &V, const example_set &E, double **I)
{
    int n=V.TotalVariables(), conse=V.Consecuente();
    int ne;

    double *px;
    double *py;
    double **pxy;
    double Inf, H;

    int nv=V.N_Antecedente();

    py = new double [V.SizeDomain(conse)];

    CalcularProbabilidadSimpleVariable(V,E,conse,V.SizeDomain(conse),py);

    cout << "PY" << endl;
    for (int i = 0; i < V.SizeDomain(conse); i++)
      cout << py[i] << " ";
    cout << endl;

    int i=0, j=0;
    double I1=0, H1=0, I2=0, H2=0, aux;
    while (i<nv)
    {
        cout << i+1 <<"/" << nv << endl;
        if (V.IsActiva(j) && V.IsAntecedente(j))
        {
            px = new double [V.SizeDomain(j)];
            pxy = new double*[V.SizeDomain(j)];
            for (int k=0; k<V.SizeDomain(j); k++)
                pxy[k]=new double [V.SizeDomain(conse)];
            CalcularProbabilidadSimpleVariable(V,E,j,V.SizeDomain(j),px);
            CalcularProbabilidadConjunta2Variables(V,E,j,conse,V.SizeDomain(j),V.SizeDomain(conse),pxy);
            I1=0;
            H1=0;
            for (int k=0; k<V.SizeDomain(conse); k++)
            {
                I2=0;
                H2=0;
                for (int q=0; q<V.SizeDomain(j); q++)
                {
                    if (pxy[q][k]==0)
                        aux=0;
                    else
                        //aux=(pxy[q][k] * log(pxy[q][k]/(px[q]*py[k])));
                        aux=pxy[q][k] * log((px[q]*py[k])/pxy[q][k]);

                    I2=I2-aux;
                    I1=I1-aux;

                    if (pxy[q][k]==0)
                        aux=0;
                    else
                        aux=pxy[q][k]*log(pxy[q][k]);

                    H2=H2-aux;
                    H1=H1-aux;
                }

                if  (H2==0)
                    I[i][k+1]=0;
                else
                    I[i][k+1]=(I2/H2);
            }
            if (H1==0)
                I[i][0]=0;
            else
                I[i][0]=(I1/H1);

            V.PrintVar(j);
            cout << ": " << I1 << "/" << H1 << ": " << I[i][0] << endl;
            for (int k=0; k<V.SizeDomain(conse); k++)
            {
                cout << "\tEtiqueta "<< k <<": " << I[i][k+1] << endl;
            }

            i++;

            delete [] px;
            for (int k=0; k<V.SizeDomain(j); k++)
                delete [] pxy[k];
            delete [] pxy;
        }
	else {
		cout << "POS: " << i << " " << j << endl;
	}
        j++;
    }

    delete [] py;


}



string extract_noun(string file)
{
    int n=file.length(), i=n-1;

    while (i>0 && file[i]!='/')
        i--;

    return file.substr(i+1,n-i-1);
}


bool Existe(string nom_fich)
{
    fstream fich(nom_fich.c_str(),ios::in);
    bool estado;

    if (!fich)
        estado = false;
    else
        estado = true;

    fich.close();

    return estado;
}



double CalcularInformacion (int tamaX, int tamaY, double *px, double *py, double **pxy)
{

    int i=0, j=0;
    double I, I1=0, H1=0, I2=0, H2=0, aux;



    for (int k=0; k<tamaY; k++)
    {
        I2=0;
        H2=0;
        for (int q=0; q<tamaX; q++)
        {
            if (pxy[q][k]==0)
                aux=0;
            else
                aux=pxy[q][k] * log((px[q]*py[k])/pxy[q][k]);

            /*cout << "I = " << I1 <<"\tH = " << H1 << endl;
            cout << "aux = " << pxy[q][k] << " log ( (" << px[q] << " * " << py[k]<<") / " << pxy[q][k] <<") = " << aux << endl;*/

            I2=I2-aux;
            I1=I1-aux;

            if (pxy[q][k]==0)
                aux=0;
            else
                aux=pxy[q][k]*log(pxy[q][k]);

            H2=H2-aux;
            H1=H1-aux;
        }

        /*if  (H2==0)
          I[i][k+1]=0;
        else
          I[i][k+1]=(I2/H2);*/
    }

    if (H1==0)
        I=0;
    else
        I=(I1/H1);

    return I;
}

double CalcularInformacion (int tamaX, int tamaY, double *px, double *py, double **pxy, double *In)
{

    int i=0, j=0;
    double I, I1=0, H1=0, I2=0, H2=0, aux;



    for (int k=0; k<tamaY; k++)
    {
        I2=0;
        H2=0;
        for (int q=0; q<tamaX; q++)
        {
            if (pxy[q][k]==0)
                aux=0;
            else
                aux=pxy[q][k] * log((px[q]*py[k])/pxy[q][k]);


            /*cout << "I = " << I1 <<"\tH = " << H1 << endl;
            cout << "aux = " << pxy[q][k] << " log ( (" << px[q] << " * " << py[k]<<") / " << pxy[q][k] <<") = " << aux << endl;*/

            if (pxy[q][k]>px[q] or pxy[q][k]>py[k]){
              cout << "La conjunta es mayor que las parciales px = " << px[q] << " py= " << py[k] <<" pxy = " << pxy[q][k] << endl;
              exit(1);
            }

            I2=I2-aux;
            I1=I1-aux;

            if (pxy[q][k]==0)
                aux=0;
            else
                aux=pxy[q][k]*log(pxy[q][k]);

            H2=H2-aux;
            H1=H1-aux;
        }

        if  (H2==0)
            In[k]=0;
        else
            In[k]=(I2/H2);
    }

    if (H1==0)
        I=0;
    else
        I=(I1/H1);

    return I;
}


double CalcularInformacion_pinta (int tamaX, int tamaY, double *px, double *py, double **pxy, double *In)
{

    int i=0, j=0;
    double I, I1=0, H1=0, I2=0, H2=0, aux;



    for (int k=0; k<tamaY; k++)
    {
        I2=0;
        H2=0;
        for (int q=0; q<tamaX; q++)
        {
            if (pxy[q][k]==0)
                aux=0;
            else
                aux=pxy[q][k] * log((px[q]*py[k])/pxy[q][k]);


            /*cout << "I = " << I1 <<"\tH = " << H1 << endl;
            cout << "aux = " << pxy[q][k] << " log ( (" << px[q] << " * " << py[k]<<") / " << pxy[q][k] <<") = " << aux << endl;*/

            if (pxy[q][k]>px[q] or pxy[q][k]>py[k]){
              cout << "La conjunta es mayor que las parciales px = " << px[q] << " py= " << py[k] <<" pxy = " << pxy[q][k] << endl;
              exit(1);
            }

            I2=I2-aux;
            I1=I1-aux;

            if (pxy[q][k]==0)
                aux=0;
            else
                aux=pxy[q][k]*log(pxy[q][k]);

            H2=H2-aux;
            H1=H1-aux;
        }

        if  (H2==0)
            In[k]=0;
        else
            In[k]=(I2/H2);
    }

    if (H1==0)
        I=0;
    else
        I=(I1/H1);

    return I;
}



void CalcularEntropia (int tamaX, int tamaY, double *px, double *py, double **pxy, double *entropia)
{

    int i=0, j=0;
    double I, I1=0, H1=0, I2=0, H2=0, aux, Iy, max=0;



    for (int k=0; k<tamaY; k++)
    {
        I2=0;
        H2=0;
        for (int q=0; q<tamaX; q++)
        {
            if (pxy[q][k]==0)
                aux=0;
            else
                aux=pxy[q][k] * log((px[q]*py[k])/pxy[q][k]);

            /*cout << "I = " << I1 <<"\tH = " << H1 << endl;
            cout << "aux = " << pxy[q][k] << " log ( (" << px[q] << " * " << py[k]<<") / " << pxy[q][k] <<") = " << aux << endl;*/

            I2=I2-aux;
            I1=I1-aux;

            if (pxy[q][k]==0)
                aux=0;
            else
                aux=(pxy[q][k]/py[k])*log(pxy[q][k]/py[k]);

            H2=H2-aux;
            H1=H1-aux;
        }

        entropia[k]=H2;
    }

    /*if (H1==0)
      I=0;
    else
      I=(I1/H1);*/

    //return max;
}




void Max_Min(const example_set &E, const VectorVar &V, int v1, int v2, int op, double &min, double &max)
{
    double result;

    if (op==-1 or op==7 or op==8)
    {
        min=E.Data(0,v1);
        max = min;
        for (int i=1; i<E.N_Examples(); i++)
        {
            result = E.Data(i,v1);
            if (result!=-999999999){
                if (min==-999999999){
                  min=result;
                  max=result;
                }
                else if (result > max)
                        max = result;
                     else if (result < min)
                        min = result;
            }
            if (result > max)
                max = result;
            else if (result < min)
                min = result;
        }

    }
    else if (op<2)
    {
        min = 0;
        max = 1;
    }
    else
    {
        /*int ss=0;
        while (ss<E.N_Examples() and Operacion(E.Data(ss,v1),E.Data(ss,v2),op)==-999999999)
          ss++;*/

        min = Operacion(E.Data(0,v1),E.Data(0,v2),op);
        max = min;
        for (int i=1; i<E.N_Examples(); i++)
        {
            double dato1, dato2;
            dato1=E.Data(i,v1);
            dato2=E.Data(i,v2);
            result = Operacion(dato1,dato2,op);
            if (result!=-999999999){
                if (min==-999999999){
                  min=result;
                  max=result;
                }
                else if (result > max)
                        max = result;
                     else if (result < min)
                        min = result;
            }
        }
    }

    //cout << "\t\t\t[" << min << "," << max <<"]\n";

}



domain_t Crear_Dominio(const VectorVar &V, example_set &E, int v1, int v2, int op)
{
    domain_t atrib;
    int n_etiquetas;
    string nombre[2]= {"Falso","Verdad"};
    double inf, sup, a[2], b[2], c[2], d[2];
    // lista de op es (0) <=, (1) >=, (2) +, (3) -, (4) dif_cuadratica, (5) *

    switch (op)
    {
    case 19: // ==
    case 20: // <=
    case 21:
    case 22:
    case 1: // >=
        n_etiquetas=2;
        inf=0;
        sup=1;
        a[0]=b[0]=c[0]=d[0]=0;
        a[1]=b[1]=c[1]=d[1]=1;
        atrib.Asigna(n_etiquetas, inf, sup, a, b, c, d, nombre);
        return atrib;
        break;
    case 2: // +
        n_etiquetas = 5;
        Max_Min(E,V,v1,v2,op,inf,sup);
        /* inf = V.Inf_Range(v1)+V.Inf_Range(v2);
            sup = V.Sup_Range(v1)+V.Sup_Range(v2);*/
        break;
    case 4: // -
        n_etiquetas = 5;
        Max_Min(E,V,v1,v2,op,inf,sup);
        break;
    case 5: // Div
        n_etiquetas = 5;
        /*inf = fabs(V.Inf_Range(v1))-V.Sup_Range(v2);
        sup = V.Sup_Range(v1)-V.Inf_Range(v2);*/
        Max_Min(E,V,v1,v2,op,inf,sup);
        break;
        /*case 4: // dif_cuadratica
                n_etiquetas = 5;
                inf =0;
                sup = pow(V.Inf_Range(v1)-V.Inf_Range(v2),2);
                if (sup < pow(V.Inf_Range(v1)-V.Sup_Range(v2),2))
                  sup = pow(V.Inf_Range(v1)-V.Sup_Range(v2),2);

                if (sup < pow(V.Sup_Range(v1)-V.Inf_Range(v2),2))
                  sup = pow(V.Sup_Range(v1)-V.Inf_Range(v2),2);

                if (sup < pow(V.Sup_Range(v1)-V.Sup_Range(v2),2))
                  sup = pow(V.Sup_Range(v1)-V.Sup_Range(v2),2);
                break;*/
    case 3: // *
        n_etiquetas = 5;
        Max_Min(E,V,v1,v2,op,inf,sup);
        /*inf = V.Inf_Range(v1)*V.Inf_Range(v2);
        sup = inf;

        if (sup < (V.Inf_Range(v1)*V.Sup_Range(v2)) )
          sup = V.Inf_Range(v1)*V.Sup_Range(v2);

        if (sup < (V.Sup_Range(v1)*V.Inf_Range(v2)) )
          sup = V.Sup_Range(v1)*V.Inf_Range(v2);

        if (sup < (V.Sup_Range(v1)*V.Sup_Range(v2)) )
          sup = V.Sup_Range(v1)*V.Sup_Range(v2);

        if (inf > (V.Inf_Range(v1)*V.Sup_Range(v2)))
               inf = (V.Inf_Range(v1)*V.Sup_Range(v2));

        if (inf > (V.Sup_Range(v1)*V.Inf_Range(v2)))
               inf = (V.Sup_Range(v1)*V.Inf_Range(v2));

        if (inf > (V.Sup_Range(v1)*V.Sup_Range(v2)))
               inf = (V.Sup_Range(v1)*V.Sup_Range(v2));

        //cout <<"Inf: " << inf << ", Sup: " << sup << endl;*/

        break;
    case 6: // Normalizacion
        n_etiquetas = 5;
        //Max_Min(E,v1,v1,-1,inf,sup);
        inf=0;
        sup=1;
        break;
    case 7: // sin
        n_etiquetas = 5;
        inf=-1;
        sup=1;
        break;
    case 8: // cos
        n_etiquetas = 5;
        inf=-1;
        sup=1;
        break;
    case 9: //log
        n_etiquetas = 5;
        Max_Min(E,V,v1,v1,9,inf,sup);
        break;
    case 10: //x^2
        n_etiquetas = 5;
        Max_Min(E,V,v1,v1,10,inf,sup);
        break;
    case 11: //sqrt
        n_etiquetas = 5;
        Max_Min(E,V,v1,v1,11,inf,sup);
        break;
    case 12: //1/
        n_etiquetas = 5;
        Max_Min(E,V,v1,v1,12,inf,sup);
        break;
    case 13: // Normalizacion
        n_etiquetas = 5;
        //Max_Min(E,v1,v1,-1,inf,sup);
        inf=-1;
        sup=1;
        break;

    case 15: //Particiones
        n_etiquetas = 5;
        Max_Min(E,V,v1,v1,-1,inf,sup);
        break;
    }
    if (inf>sup)
    {
        cout << "ERROR!!! en el rango de: " << v1 << " " << v2 << " con la operacion " << Simbolo_Operador(op) << endl;
        exit(1);
    }

    if (op==4 and (inf<0 and sup>0))
      atrib.Asigna_Op_menos(n_etiquetas,inf,sup,false,false);
    else
      atrib.Asigna(n_etiquetas,inf,sup,false,false);

    return atrib;
}

domain_t Crear_Dominio(const VectorVar &V, example_set &E, int v1, int v2, int op, int n_labels)
{
    domain_t atrib;
    int n_etiquetas;
    string nombre[2]= {"Falso","Verdad"};
    double inf, sup, a[2], b[2], c[2], d[2];
    // lista de op es (0) <=, (1) >=, (2) +, (3) -, (4) dif_cuadratica, (5) *

    switch (op)
    {
    case 19: // ==
    case 20: // <=
    case 21:
    case 22:
    case 1: // >=
        n_etiquetas=2;
        inf=0;
        sup=1;
        a[0]=b[0]=c[0]=d[0]=0;
        a[1]=b[1]=c[1]=d[1]=1;
        atrib.Asigna(n_etiquetas, inf, sup, a, b, c, d, nombre);
        return atrib;
        break;
    case 2: // +
        n_etiquetas = n_labels;
        Max_Min(E,V,v1,v2,op,inf,sup);
        /* inf = V.Inf_Range(v1)+V.Inf_Range(v2);
            sup = V.Sup_Range(v1)+V.Sup_Range(v2);*/
        break;
    case 4:
        n_etiquetas = n_labels;
        Max_Min(E,V,v1,v2,op,inf,sup);
        break;
    case 5: // -
        n_etiquetas = n_labels;
        /*inf = fabs(V.Inf_Range(v1))-V.Sup_Range(v2);
        sup = V.Sup_Range(v1)-V.Inf_Range(v2);*/
        Max_Min(E,V,v1,v2,op,inf,sup);
        break;
        /*case 4: // dif_cuadratica
                n_etiquetas = 5;
                inf =0;
                sup = pow(V.Inf_Range(v1)-V.Inf_Range(v2),2);
                if (sup < pow(V.Inf_Range(v1)-V.Sup_Range(v2),2))
                  sup = pow(V.Inf_Range(v1)-V.Sup_Range(v2),2);

                if (sup < pow(V.Sup_Range(v1)-V.Inf_Range(v2),2))
                  sup = pow(V.Sup_Range(v1)-V.Inf_Range(v2),2);

                if (sup < pow(V.Sup_Range(v1)-V.Sup_Range(v2),2))
                  sup = pow(V.Sup_Range(v1)-V.Sup_Range(v2),2);
                break;*/
    case 3: // *
        n_etiquetas = n_labels;
        Max_Min(E,V,v1,v2,op,inf,sup);
        /*inf = V.Inf_Range(v1)*V.Inf_Range(v2);
        sup = inf;

        if (sup < (V.Inf_Range(v1)*V.Sup_Range(v2)) )
          sup = V.Inf_Range(v1)*V.Sup_Range(v2);

        if (sup < (V.Sup_Range(v1)*V.Inf_Range(v2)) )
          sup = V.Sup_Range(v1)*V.Inf_Range(v2);

        if (sup < (V.Sup_Range(v1)*V.Sup_Range(v2)) )
          sup = V.Sup_Range(v1)*V.Sup_Range(v2);

        if (inf > (V.Inf_Range(v1)*V.Sup_Range(v2)))
               inf = (V.Inf_Range(v1)*V.Sup_Range(v2));

        if (inf > (V.Sup_Range(v1)*V.Inf_Range(v2)))
               inf = (V.Sup_Range(v1)*V.Inf_Range(v2));

        if (inf > (V.Sup_Range(v1)*V.Sup_Range(v2)))
               inf = (V.Sup_Range(v1)*V.Sup_Range(v2));

        //cout <<"Inf: " << inf << ", Sup: " << sup << endl;*/

        break;
    case 6: // Normalizacion
        n_etiquetas = n_labels;
        Max_Min(E,V,v1,v1,-1,inf,sup);
        inf=0;
        sup=1;
        break;
    case 7: // sin
        Max_Min(E,V,v1,v1,op,inf,sup);
        n_etiquetas = n_labels;
        break;
    case 8: // cos
        Max_Min(E,V,v1,v1,op,inf,sup);
        n_etiquetas = n_labels;
        break;
    case 9: //log
        n_etiquetas = n_labels;
        Max_Min(E,V,v1,v1,9,inf,sup);
        break;
    case 10: //X^2
        n_etiquetas = n_labels;
        Max_Min(E,V,v1,v1,10,inf,sup);
        break;
    case 11: //sqrt
        n_etiquetas = n_labels;
        Max_Min(E,V,v1,v1,11,inf,sup);
        break;
    case 12: // 1/x
        n_etiquetas = n_labels;
        Max_Min(E,V,v1,v1,12,inf,sup);
        break;
    case 13: // Normalizacion Extendida
        n_etiquetas = n_labels;
        Max_Min(E,V,v1,v1,-1,inf,sup);
        inf=-1;
        sup=1;
        break;

    case 15: //Particiones
        n_etiquetas = n_labels;
        Max_Min(E,V,v1,v1,-1,inf,sup);
        break;
    }
    if (inf>sup)
    {
        cout << "ERROR!!! en el rango de: " << v1 << " " << v2 << " con la operacion " << Simbolo_Operador(op) << endl;
        exit(1);
    }

    if (op==4 and (inf<0 and sup>0))
      atrib.Asigna_Op_menos(n_etiquetas,inf,sup,false,false);
    else
      atrib.Asigna(n_etiquetas,inf,sup,false,false);
    return atrib;
}



domain_t Crear_Dominio(const VectorVar &V, example_set &E, int v1, int v2, int op, int n_labels, int desp)
{
    domain_t atrib;
    int n_etiquetas;
    string nombre[2]= {"Falso","Verdad"};
    double inf, sup, a[2], b[2], c[2], d[2];
    // lista de op es (0) <=, (1) >=, (2) +, (3) -, (4) dif_cuadratica, (5) *

    switch (op)
    {
    case 19: // ==
    case 20: // <=
    case 21:
    case 22:
    case 1: // >=
        n_etiquetas=2;
        inf=0;
        sup=1;
        a[0]=b[0]=c[0]=d[0]=0;
        a[1]=b[1]=c[1]=d[1]=1;
        atrib.Asigna(n_etiquetas, inf, sup, a, b, c, d, nombre);
        return atrib;
        break;
    case 2: // +
        n_etiquetas = n_labels;
        Max_Min(E,V,v1,v2,op,inf,sup);
        /* inf = V.Inf_Range(v1)+V.Inf_Range(v2);
            sup = V.Sup_Range(v1)+V.Sup_Range(v2);*/
        break;
    case 4:
        n_etiquetas = n_labels;
        Max_Min(E,V,v1,v2,op,inf,sup);
        break;
    case 5: // -
        n_etiquetas = n_labels;
        /*inf = fabs(V.Inf_Range(v1))-V.Sup_Range(v2);
        sup = V.Sup_Range(v1)-V.Inf_Range(v2);*/
        Max_Min(E,V,v1,v2,op,inf,sup);
        break;
        /*case 4: // dif_cuadratica
                n_etiquetas = 5;
                inf =0;
                sup = pow(V.Inf_Range(v1)-V.Inf_Range(v2),2);
                if (sup < pow(V.Inf_Range(v1)-V.Sup_Range(v2),2))
                  sup = pow(V.Inf_Range(v1)-V.Sup_Range(v2),2);

                if (sup < pow(V.Sup_Range(v1)-V.Inf_Range(v2),2))
                  sup = pow(V.Sup_Range(v1)-V.Inf_Range(v2),2);

                if (sup < pow(V.Sup_Range(v1)-V.Sup_Range(v2),2))
                  sup = pow(V.Sup_Range(v1)-V.Sup_Range(v2),2);
                break;*/
    case 3: // *
        n_etiquetas = n_labels;
        Max_Min(E,V,v1,v2,op,inf,sup);
        /*inf = V.Inf_Range(v1)*V.Inf_Range(v2);
        sup = inf;

        if (sup < (V.Inf_Range(v1)*V.Sup_Range(v2)) )
          sup = V.Inf_Range(v1)*V.Sup_Range(v2);

        if (sup < (V.Sup_Range(v1)*V.Inf_Range(v2)) )
          sup = V.Sup_Range(v1)*V.Inf_Range(v2);

        if (sup < (V.Sup_Range(v1)*V.Sup_Range(v2)) )
          sup = V.Sup_Range(v1)*V.Sup_Range(v2);

        if (inf > (V.Inf_Range(v1)*V.Sup_Range(v2)))
               inf = (V.Inf_Range(v1)*V.Sup_Range(v2));

        if (inf > (V.Sup_Range(v1)*V.Inf_Range(v2)))
               inf = (V.Sup_Range(v1)*V.Inf_Range(v2));

        if (inf > (V.Sup_Range(v1)*V.Sup_Range(v2)))
               inf = (V.Sup_Range(v1)*V.Sup_Range(v2));

        //cout <<"Inf: " << inf << ", Sup: " << sup << endl;*/

        break;
    case 6: // Normalizacion
        n_etiquetas = n_labels;
        //Max_Min(E,v1,v1,-1,inf,sup);
        inf=0;
        sup=1;
        break;
    case 7: // Delta+
        n_etiquetas = n_labels;
        Max_Min(E,V,v1,v1,op,inf,sup);
        //cout <<"sin [" << inf << "," << sup <<"]\n";
        break;
    case 8: // Delta-
        n_etiquetas = n_labels;
        Max_Min(E,V,v1,v1,op,inf,sup);
        //cout <<"cos [" << inf << "," << sup <<"]\n";
        break;
    case 9: //log
        n_etiquetas = n_labels;
        Max_Min(E,V,v1,v1,9,inf,sup);
        break;
    case 10: //X^2
        n_etiquetas = n_labels;
        Max_Min(E,V,v1,v1,10,inf,sup);
        break;
    case 11: //sqrt
        n_etiquetas = n_labels;
        Max_Min(E,V,v1,v1,11,inf,sup);
        break;
    case 12: // 1/x
        n_etiquetas = n_labels;
        Max_Min(E,V,v1,v1,12,inf,sup);
        break;
    case 13: // Normalizacion Extendida
        n_etiquetas = n_labels;
        //Max_Min(E,v1,v1,-1,inf,sup);
        inf=-1;
        sup=1;
        break;
    case 15: //Particiones
        n_etiquetas = n_labels;
        Max_Min(E,V,v1,v1,-1,inf,sup);
        break;
    }
    if (inf>sup)
    {
        cout << "ERROR!!! en el rango de: " << v1 << " " << v2 << " con la operacion " << Simbolo_Operador(op) << endl;
        exit(1);
    }

    /*if (op==4 and (inf<0 and sup>0))
      atrib.Asigna_Op_menos(n_etiquetas,inf,sup,false,false);
    else*/

   if (op!=7 and op!=8)
      atrib.Asigna(n_etiquetas,inf,sup,false,false,desp);
   else if (op==7){
          inf=inf-((sup-inf)/(4*n_etiquetas));
          atrib.Asigna(n_etiquetas,inf,sup,false,false,desp);
        }
        else {
          sup=sup+((sup-inf)/(4*n_etiquetas));
          atrib.Asigna(n_etiquetas,inf,sup,false,false,desp);
        }
    return atrib;
}




string Simbolo_Operador(int op)
{
    string cadena;
    switch (op)
    {
    case 19:
	cadena = " == ";
	break;
    case 20:
        cadena = " OR ";
        break;
    case 21:
        cadena = " XOR ";
        break;
    case 22:
        cadena = " AND ";
        break;
    case 1:
        cadena = ">=";
        break;
    case 2:
        cadena = "+";
        break;
    case 5:
        cadena = "/";
        break;
    case 4:
        cadena = "-";
        break;
    case 3:
        cadena = "*";
        break;
    case 6:
        cadena = "nor";
        break;
    case 7:
        cadena = "Delta+";
        break;
    case 8:
        cadena = "Delta-";
        break;
    case 9:
        cadena = "log";
        break;
    case 10:
        cadena = "pow";
        break;
    case 11:
        cadena = "sqrt";
        break;
    case 12:
        cadena = "inv";
        break;
    case 13:
        cadena = "nor_ext";
        break;
    case 15:
        cadena = "Domain";
        break;
    }
    return cadena;
}



void InsertarOrdenado(double info[], int orden[], int tama)
{
    int i=tama-1, aux;
    orden[tama]=tama;
    while (i>=0 && info[orden[i]]<info[orden[i+1]])
    {
        aux=orden[i];
        orden[i]=orden[i+1];
        orden[i+1]=aux;
        i--;
    }
}


// Calculo de la ganancia de informacion de nuevos atributos creados a partir de operadores
void MedidasDeInformacion_operadores(VectorVar &V, example_set &E, VectorVar &F)
{
    struct nova_var
    {
        int v1;
        int v2;
        int op;
    };

    const char *intlabel[6]= {"2","3","4","5","6","7"};
    const char *intsensibilidad[4]= {"Uni","Izq","Dch","Cnt"};
    int n=V.TotalVariables(), conse=V.Consecuente();
    int ne;
    bool anadir;
    int ciclo;

    double *px;
    double *pz;
    double *py;
    double *pf;
    double **pxy;
    double **pzy;
    double **pfy;
    double Inf, H;
    double *Information, **entropia;
    int *orden;
    int **orden_clase;

    double measure[100];

    double minimo[V.SizeDomain(conse)],/* entropia[V.SizeDomain(conse)],*/ max_clase[V.SizeDomain(conse)],
    nuevo_max[V.SizeDomain(conse)];

    nova_var a_introducir[V.SizeDomain(conse)];
    string NomVar1, NomVar2, NomVar3;
    domain_t atrib;

    for (int i=0; i<V.SizeDomain(conse); i++)
    {
        max_clase[i]=0;
        minimo[i]=100000.0;
        a_introducir[i].v1=-1;
        a_introducir[i].v2=-1;
        a_introducir[i].op=-1;
    }

    int nv=V.N_Antecedente();
    Information = new double[nv];
    entropia = new double*[V.SizeDomain(conse)];
    orden = new int[nv];
    orden_clase = new int*[V.SizeDomain(conse)];
    for (int i=0; i<V.SizeDomain(conse); i++)
    {
        orden_clase[i]= new int[nv];
        entropia[i]= new double [nv];
    }
    int nmax=0;

    py = new double [V.SizeDomain(conse)];

    CalcularProbabilidadSimpleVariable(V,E,conse,V.SizeDomain(conse),py);

    double max=0, min=10000;
    for (int x=0; x<nv; x++)
    {
        px = new double [V.SizeDomain(x)];
        CalcularProbabilidadSimpleVariable(V,E,x,V.SizeDomain(x),px);
        pxy = new double*[V.SizeDomain(x)];
        for (int k=0; k<V.SizeDomain(x); k++)
            pxy[k]=new double [V.SizeDomain(conse)];
        CalcularProbabilidadConjunta2Variables(V,E,x,conse,V.SizeDomain(x),V.SizeDomain(conse),pxy);

        Information[x] = CalcularInformacion(V.SizeDomain(x), V.SizeDomain(conse), px, py, pxy, measure);
        InsertarOrdenado (Information, orden, x);
        //CalcularEntropia(V.SizeDomain(x), V.SizeDomain(conse), px, py, pxy, entropia);

        for (int k=0; k<V.SizeDomain(x); k++)
            delete [] pxy[k];
        delete [] pxy;
        delete [] px;

        //for (int i=0; i<V.SizeDomain(conse); i++)
        if (Information[x]>max)
            max=Information[x]; //+0.0001;

        NomVar1 = V.SPrintVar(x);

        cout << "Informacion (" << NomVar1 << ")= " << Information[x]<<"     ";
        for (int i=0; i<V.SizeDomain(conse); i++)
        {
            entropia[i][x]=measure[i];
            InsertarOrdenado (entropia[i], orden_clase[i],x);
            cout << entropia[i][x] << " | ";
            if (entropia[i][x]>max_clase[i])
            {
                max_clase[i]=entropia[i][x];
                nuevo_max[i]=max_clase[i];
            }
        }
        cout << endl;
    }

    cout << "--------------------------------------------------------------------------------------------" << endl;
    for (int i=0; i<V.SizeDomain(conse); i++)
        cout << max_clase[i] << " | ";
    cout << endl;
    cout << "--------------------------------------------------------------------------------------------" << endl;

    // Trabajar solo con las 100 variables mas informadas
    cout << "--------------------------------------------------------------------------------------------" << endl;
    for (int i=0; i<nv && i<100; i++)
        cout << Information[orden[i]] << endl;
    cout << endl;
    cout << "--------------------------------------------------------------------------------------------" << endl;


    for (int x=0; x<nv-1; x++)
    {
      if (V.Get_Unit(x)!=0) {
        NomVar1= V.SPrintVar(x);
        cout <<"----------------------------------------------------  " << NomVar1 << "--- "<< V.Get_Unit(x) << " ----------------------------\n";


        for (int z=x+1; z<nv; z++)
        {
            if (/*!V.IsDiscrete(z) and*/ V.Get_Unit(z)!=0)
                for (int op=19; op<=19 /*6 */; op++){
                    if  ( (V.Get_Unit(z)==V.Get_Unit(x)) ){
                       atrib = Crear_Dominio(V,E,x,z,op);
                       pf = new double [atrib.Size()];
                       CalcularProbabilidadSimpleVariable_op(V,E,x,z,atrib,op,pf);
                       pfy = new double*[atrib.Size()];
                       for (int k=0; k<atrib.Size(); k++)
                         pfy[k]=new double [V.SizeDomain(conse)];

                       CalcularProbabilidadConjunta2Variables_op(V,E,x,z,atrib,conse,V.SizeDomain(conse),op, pfy);

                       Inf = CalcularInformacion(atrib.Size(), V.SizeDomain(conse), pf, py, pfy, measure);

                       bool mejora=false;
                       for (int i=0; i<V.SizeDomain(conse) /*&& !mejora*/; i++){
                         if (measure[i]>(0.2*max_clase[i])){
                            mejora = true;
                            cout << "mejora para la clase " << i << " (de " << max_clase[i] << " a " << measure[i] <<")\n";
                            //cout << "Nuevo valor: " << entropia[i] << "   Anterior Valor: " << max_clase[i] << endl;
                         }
                       }



                       if (mejora && conse<999){
                          nmax++;

                          NomVar1 = V.SPrintVar(x);
                          NomVar2 = V.SPrintVar(z);
                          NomVar3 = "("+NomVar1+Simbolo_Operador(op)+NomVar2+")";
                          cout << NomVar3 << "-> "<< Inf << "( " << measure[0];

                          for (int i=1; i<V.SizeDomain(conse) /*&& !mejora*/; i++)
                            cout << ", " << measure[i];
                          cout <<" )"<< endl;


                          variable_t nueva_variable(NomVar3,true,true,atrib);

                          nueva_variable.Pinta();
                          F.Add_Variable_Lista(nueva_variable,op,x,z,nueva_variable.Inf_Range(),nueva_variable.Sup_Range());

                          for (int i=0; i<V.SizeDomain(conse); i++)
                            if (measure[i]>max_clase[i])
                              nuevo_max[i]=measure[i];

                      }
                      for (int k=0; k<atrib.Size(); k++)
                         delete [] pfy[k];
                      delete [] pfy;
                      delete [] pf;


                      if (op==5) {
                       atrib = Crear_Dominio(V,E,z,x,op,9,0);
                       pf = new double [atrib.Size()];
                       CalcularProbabilidadSimpleVariable_op(V,E,z,x,atrib,op,pf);
                       pfy = new double*[atrib.Size()];
                       for (int k=0; k<atrib.Size(); k++)
                         pfy[k]=new double [V.SizeDomain(conse)];

                       CalcularProbabilidadConjunta2Variables_op(V,E,z,x,atrib,conse,V.SizeDomain(conse),op, pfy);

                       Inf = CalcularInformacion(atrib.Size(), V.SizeDomain(conse), pf, py, pfy, measure);

                       bool mejora=false;
                       for (int i=0; i<V.SizeDomain(conse) /*&& !mejora*/; i++){
                         if (measure[i]>(0.9*max_clase[i])){
                            mejora = true;
                            cout << "mejora para la clase " << i << " (de " << max_clase[i] << " a " << measure[i] <<")\n";
                            //cout << "Nuevo valor: " << entropia[i] << "   Anterior Valor: " << max_clase[i] << endl;
                         }
                       }




                       if (mejora && conse<999){
                          nmax++;

                          NomVar1 = V.SPrintVar(z);
                          NomVar2 = V.SPrintVar(x);
                          NomVar3 = "("+NomVar1+Simbolo_Operador(op)+NomVar2+")";
                          cout << NomVar3 << "-> "<< Inf << "( " << measure[0];

                          for (int i=1; i<V.SizeDomain(conse) /*&& !mejora*/; i++)
                            cout << ", " << measure[i];
                          cout <<" )"<< endl;

                          variable_t nueva_variable(NomVar3,true,true,atrib);
                          F.Add_Variable_Lista(nueva_variable,op,z,x,nueva_variable.Inf_Range(),nueva_variable.Sup_Range());

                          for (int i=0; i<V.SizeDomain(conse); i++)
                            if (measure[i]>max_clase[i])
                              nuevo_max[i]=measure[i];

                        }
                        for (int k=0; k<atrib.Size(); k++)
                          delete [] pfy[k];
                        delete [] pfy;
                        delete [] pf;

                      } // op==5



                }

            }

        }
      }
    }


    cout << "Numero de atributos prometedores= " << nmax << endl;
    //cin.get();
    for (int i=0; i<V.SizeDomain(conse); i++)
        delete [] orden_clase[i];
    delete [] orden_clase;

    delete [] py;
    delete [] Information;
    delete [] orden;
}


// Determina la mejor particion del dominio de una funcion

domain_t Mejor_Dominio_Operacion(VectorVar &V, example_set &E, int x, int z,int op, double *py, int conse, int &sensi){
  domain_t atrib;
  double *pf, **pfy;
  double MaxInf=0, Inf;
  int mejor_etiqs=2, mejor_sensib=0;
  bool primera_vez=true;

  double *measure = new double[V.SizeDomain(conse)];
  double *mejor_measure = new double[V.SizeDomain(conse)];
  double *basic_measure = new double[V.SizeDomain(conse)];

  for (int i=0; i<V.SizeDomain(conse); i++)
    mejor_measure[i]=0;

  for (int etiq=5; etiq<=5; etiq++){
    for (int sensibilidad=0; sensibilidad<1; sensibilidad++){
       atrib = Crear_Dominio(V,E,x,z,op,etiq,sensibilidad);
       pf = new double [atrib.Size()];
       CalcularProbabilidadSimpleVariable_op(V,E,x,z,atrib,op,pf);
       pfy = new double*[atrib.Size()];
       for (int k=0; k<atrib.Size(); k++)
         pfy[k]=new double [V.SizeDomain(conse)];

       CalcularProbabilidadConjunta2Variables_op(V,E,x,z,atrib,conse,V.SizeDomain(conse),op, pfy);

       Inf = CalcularInformacion(atrib.Size(), V.SizeDomain(conse), pf, py, pfy, measure);

       /*if (primera_vez){
          primera_vez=false;
          for (int i=0; i<V.SizeDomain(conse); i++){
            mejor_measure[i]=1;
            basic_measure[i]=measure[i];
            Inf=0.1;
          }
       }
       else {
          int conteo_basic=0, conteo_mejor=0;
          for (int i=0; i<V.SizeDomain(conse); i++){
            if (measure[i]>=basic_measure[i])
              conteo_basic++;
              if ( (measure[i]+1)/(basic_measure[i]+1) > mejor_measure[i] ){
                 mejor_measure[i]= (measure[i]+1)/(basic_measure[i]+1);
                 conteo_mejor++;
              }
          }
          Inf=conteo_basic+(conteo_mejor*0.1);

       }*/


       if (Inf>MaxInf){
         MaxInf=Inf;
         mejor_etiqs=etiq;
         mejor_sensib=sensibilidad;
       }

       for (int k=0; k<atrib.Size(); k++)
         delete [] pfy[k];
       delete [] pfy;
       delete [] pf;
    }
  }

  atrib = Crear_Dominio(V,E,x,z,op,mejor_etiqs, mejor_sensib);
  sensi = mejor_sensib;


  delete [] measure;
  delete [] mejor_measure;
  delete [] basic_measure;

  return atrib;

}



    class nodo{
       public:
       double valor;
       int operando1;
       int operando2;
       int operacion;
       int labels;
       int sensibilidad;
       domain_t atributo;
       string nombre;
       bool operator< (const nodo &nodo2) const {
         return (valor<nodo2.valor);
       }
       nodo (){
       }
       nodo(const nodo &nodo2){
         valor=nodo2.valor;
         operando1=nodo2.operando1;
         operando2=nodo2.operando2;
         operacion=nodo2.operacion;
         labels=nodo2.labels;
         sensibilidad=nodo2.sensibilidad;
         atributo=nodo2.atributo;
         nombre=nodo2.nombre;
       }
    };




// Calculo de la ganancia de informacion de nuevos atributos creados a partir de operadores ajustando el dominio
void MedidasDeInformacion_operadores_AjustandoDominio(VectorVar &V, example_set &E, VectorVar &F)
{
    struct nova_var
    {
        int v1;
        int v2;
        int op;
    };


    priority_queue<nodo, vector<nodo>, less<nodo> > mi_cola;

    const char *intlabel[6]= {"2","3","4","5","6","7"};
    const char *intsensibilidad[4]= {"Uni","Izq","Dch","Cnt"};
    int n=V.TotalVariables(), conse=V.Consecuente();
    int ne;
    bool anadir;
    int ciclo;

    double *px;
    double *pz;
    double *py;
    double *pf;
    double **pxy;
    double **pzy;
    double **pfy;
    double Inf, H;
    double *Information, **entropia;
    int *orden;
    int **orden_clase;
    nodo candidato;

    double measure[100];

    double minimo[V.SizeDomain(conse)],/* entropia[V.SizeDomain(conse)],*/ max_clase[V.SizeDomain(conse)],
    nuevo_max[V.SizeDomain(conse)];

    nova_var a_introducir[V.SizeDomain(conse)];
    string NomVar1, NomVar2, NomVar3;
    domain_t atrib, atrib2;

    for (int i=0; i<V.SizeDomain(conse); i++)
    {
        max_clase[i]=0;
        minimo[i]=100000.0;
        a_introducir[i].v1=-1;
        a_introducir[i].v2=-1;
        a_introducir[i].op=-1;
    }

    int nv=V.N_Antecedente();
    Information = new double[nv];
    entropia = new double*[V.SizeDomain(conse)];
    orden = new int[nv];
    orden_clase = new int*[V.SizeDomain(conse)];
    for (int i=0; i<V.SizeDomain(conse); i++)
    {
        orden_clase[i]= new int[nv];
        entropia[i]= new double [nv];
    }
    int nmax=0;

    py = new double [V.SizeDomain(conse)];

    CalcularProbabilidadSimpleVariable(V,E,conse,V.SizeDomain(conse),py);

    double max=0, min=10000;
    // Variables antecedente de V
    for (int x=0; x<nv; x++)
    {
        px = new double [V.SizeDomain(x)];
        CalcularProbabilidadSimpleVariable(V,E,x,V.SizeDomain(x),px);
        pxy = new double*[V.SizeDomain(x)];
        for (int k=0; k<V.SizeDomain(x); k++)
            pxy[k]=new double [V.SizeDomain(conse)];
        CalcularProbabilidadConjunta2Variables(V,E,x,conse,V.SizeDomain(x),V.SizeDomain(conse),pxy);

        Information[x] = CalcularInformacion(V.SizeDomain(x), V.SizeDomain(conse), px, py, pxy, measure);
        InsertarOrdenado (Information, orden, x);
        //CalcularEntropia(V.SizeDomain(x), V.SizeDomain(conse), px, py, pxy, entropia);

        for (int k=0; k<V.SizeDomain(x); k++)
            delete [] pxy[k];
        delete [] pxy;
        delete [] px;

        //for (int i=0; i<V.SizeDomain(conse); i++)
        if (Information[x]>max)
            max=Information[x]; //+0.0001;

        NomVar1 = V.SPrintVar(x);

        cout << "Informacion (" << NomVar1 << ")= " << Information[x]<<"     ";
        for (int i=0; i<V.SizeDomain(conse); i++)
        {
            entropia[i][x]=measure[i];
            InsertarOrdenado (entropia[i], orden_clase[i],x);
            cout << entropia[i][x] << " | ";
            if (entropia[i][x]>max_clase[i])
            {
                max_clase[i]=entropia[i][x];
                nuevo_max[i]=max_clase[i];
            }
        }
        cout << endl;
    }

    cout << "--------------------------------------------------------------------------------------------" << endl;
    for (int i=0; i<V.SizeDomain(conse); i++)
        cout << max_clase[i] << " | ";
    cout << endl;
    cout << "--------------------------------------------------------------------------------------------" << endl;

    // Trabajar solo con las 100 variables mas informadas
    cout << "--------------------------------------------------------------------------------------------" << endl;
    for (int i=0; i<nv && i<100; i++)
        cout << Information[orden[i]] << endl;
    cout << endl;
    cout << "--------------------------------------------------------------------------------------------" << endl;

    int sensi;

    for (int x=0; x<nv; x++)
    {
      if (V.Get_Unit(x)!=0) {
        NomVar1= V.SPrintVar(x);
        cout <<"----------------------------------------------------  " << NomVar1 << "--- "<< V.Get_Unit(x) << " ----------------------------\n";

        if (!V.IsDiscrete(x)){
            atrib = V.Domain(x);
            atrib.Amplia_Delta();
            pf = new double [atrib.Size()];
            CalcularProbabilidadSimpleVariable_op(V,E,x,x,atrib,7,pf);
            pfy = new double*[atrib.Size()];
            for (int k=0; k<atrib.Size(); k++)
            pfy[k]=new double [V.SizeDomain(conse)];

            CalcularProbabilidadConjunta2Variables_op(V,E,x,x,atrib,conse,V.SizeDomain(conse),7, pfy);

            Inf = CalcularInformacion(atrib.Size(), V.SizeDomain(conse), pf, py, pfy, measure);

            NomVar1 = V.SPrintVar(x);
            NomVar2 = V.SPrintVar(x);
            NomVar3 = Simbolo_Operador(7)+"("+NomVar1+")";
            cout << NomVar3 << "-> "<< Inf << "( " << measure[0];
                          for (int i=1; i<V.SizeDomain(conse); i++)
                            cout << ", " << measure[i];
                          cout <<" )"<< endl;

            candidato.valor = Inf;
            candidato.operando1 = x;
            candidato.operando2 = x;
            candidato.operacion = 7;
            candidato.nombre = NomVar3;
            candidato.labels = atrib.Size();
            candidato.sensibilidad = sensi;
            candidato.atributo = atrib;

            mi_cola.push(candidato);
            nmax++;

            atrib = V.Domain(x);
            atrib.Reduce_Delta();
            CalcularProbabilidadSimpleVariable_op(V,E,x,x,atrib,8,pf);
            CalcularProbabilidadConjunta2Variables_op(V,E,x,x,atrib,conse,V.SizeDomain(conse),8, pfy);

            Inf = CalcularInformacion(atrib.Size(), V.SizeDomain(conse), pf, py, pfy, measure);

            NomVar1 = V.SPrintVar(x);
            NomVar2 = V.SPrintVar(x);
            NomVar3 = Simbolo_Operador(8)+"("+NomVar1+")";
            cout << NomVar3 << "-> "<< Inf << "( " << measure[0];
                          for (int i=1; i<V.SizeDomain(conse); i++)
                            cout << ", " << measure[i];
                          cout <<" )"<< endl;

            candidato.valor = Inf;
            candidato.operando1 = x;
            candidato.operando2 = x;
            candidato.operacion = 8;
            candidato.nombre = NomVar3;
            candidato.labels = atrib.Size();
            candidato.sensibilidad = sensi;
            candidato.atributo = atrib;

            mi_cola.push(candidato);
            nmax++;


            for (int k=0; k<atrib.Size(); k++)
              delete [] pfy[k];
            delete [] pfy;
            delete [] pf;

        }


        for (int z=x+1; z<nv; z++)
        {
            if ( V.Get_Unit(z)!=0)
                for (int op=1; op<5; op++){
                    if  ( ( (op==1 or op==2 or op==4) and (V.Get_Unit(z)==V.Get_Unit(x)) ) or op==3 or (op==5 )){
                       if (op!=1 )
                         atrib = Mejor_Dominio_Operacion(V,E,x,z,op, py, conse, sensi);
                       else
                         atrib = Crear_Dominio(V,E,x,z,op,2,0);

                       pf = new double [atrib.Size()];
                       CalcularProbabilidadSimpleVariable_op(V,E,x,z,atrib,op,pf);
                       pfy = new double*[atrib.Size()];
                       for (int k=0; k<atrib.Size(); k++)
                         pfy[k]=new double [V.SizeDomain(conse)];

                       CalcularProbabilidadConjunta2Variables_op(V,E,x,z,atrib,conse,V.SizeDomain(conse),op, pfy);

                       Inf = CalcularInformacion(atrib.Size(), V.SizeDomain(conse), pf, py, pfy, measure);

                       bool mejora=false;
                       double temp=0;
                       for (int i=0; i<V.SizeDomain(conse); i++){
                         if (measure[i]>((1.0*max_clase[i]))){
                            if ( (measure[i]/max_clase[i]) > temp)
                              //temp=measure[i]/max_clase[i];
                              temp=Inf;
                            mejora = true;
                            cout << "mejora para la clase " << i << " (de " << max_clase[i] << " a " << measure[i] <<")\n";
                            //cout << "Nuevo valor: " << entropia[i] << "   Anterior Valor: " << max_clase[i] << endl;
                         }
                       }



                       if (mejora && conse<999){
                          nmax+=3;

                          NomVar1 = V.SPrintVar(x);
                          NomVar2 = V.SPrintVar(z);
                          NomVar3 = "("+NomVar1+Simbolo_Operador(op)+NomVar2+")";
                          cout << NomVar3 << "-> "<< Inf << "( " << measure[0];
                          for (int i=1; i<V.SizeDomain(conse); i++)
                            cout << ", " << measure[i];
                          cout <<" )"<< endl;

                          candidato.valor = Inf;
                          candidato.operando1 = x;
                          candidato.operando2 = z;
                          candidato.operacion = op;
                          candidato.nombre = NomVar3;
                          candidato.labels = atrib.Size();
                          candidato.sensibilidad = sensi;
                          candidato.atributo = atrib;

                          mi_cola.push(candidato);

                       if (op!=1){
                          // Meter Desplazamiento del dominio de crecimiento
                          string Delta="Delta+";
                          atrib2 = atrib;
                          atrib2.Amplia_Delta();
                          CalcularProbabilidadSimpleVariable_op(V,E,x,z,atrib2,op,pf);
                          CalcularProbabilidadConjunta2Variables_op(V,E,x,z,atrib2,conse,V.SizeDomain(conse),op, pfy);
                          Inf = CalcularInformacion(atrib2.Size(), V.SizeDomain(conse), pf, py, pfy, measure);

                          NomVar1 = V.SPrintVar(x);
                          NomVar2 = V.SPrintVar(z);
                          NomVar3 = Delta+"("+NomVar1+Simbolo_Operador(op)+NomVar2+")";
                          cout << NomVar3 << "-> "<< Inf << "( " << measure[0];
                          for (int i=1; i<V.SizeDomain(conse); i++)
                            cout << ", " << measure[i];
                          cout <<" )"<< endl;

                          candidato.valor = Inf;
                          candidato.operando1 = x;
                          candidato.operando2 = z;
                          candidato.operacion = op;
                          candidato.nombre = NomVar3;
                          candidato.labels = atrib2.Size();
                          candidato.sensibilidad = sensi;
                          candidato.atributo = atrib2;

                          mi_cola.push(candidato);

                          // Meter Desplazamiento del dominio de decrecimiento
                          Delta="Delta-";
                          atrib2 = atrib;
                          atrib2.Reduce_Delta();
                          CalcularProbabilidadSimpleVariable_op(V,E,x,z,atrib2,op,pf);
                          CalcularProbabilidadConjunta2Variables_op(V,E,x,z,atrib2,conse,V.SizeDomain(conse),op, pfy);
                          Inf = CalcularInformacion(atrib2.Size(), V.SizeDomain(conse), pf, py, pfy, measure);

                          NomVar1 = V.SPrintVar(x);
                          NomVar2 = V.SPrintVar(z);
                          NomVar3 = Delta+"("+NomVar1+Simbolo_Operador(op)+NomVar2+")";
                          cout << NomVar3 << "-> "<< Inf << "( " << measure[0];
                          for (int i=1; i<V.SizeDomain(conse); i++)
                            cout << ", " << measure[i];
                          cout <<" )"<< endl;

                          candidato.valor = Inf;
                          candidato.operando1 = x;
                          candidato.operando2 = z;
                          candidato.operacion = op;
                          candidato.nombre = NomVar3;
                          candidato.labels = atrib2.Size();
                          candidato.sensibilidad = sensi;
                          candidato.atributo = atrib2;

                          mi_cola.push(candidato);
                        }


                          for (int i=0; i<V.SizeDomain(conse); i++)
                            if (measure[i]>max_clase[i])
                              nuevo_max[i]=measure[i];

                      }
                      for (int k=0; k<atrib.Size(); k++)
                         delete [] pfy[k];
                      delete [] pfy;
                      delete [] pf;





                }

            }

        }
      }
    }


    cout << "Numero de atributos prometedores= " << nmax << endl;
    //cin.get();
    for (int i=0; i<V.SizeDomain(conse); i++)
        delete [] orden_clase[i];
    delete [] orden_clase;


    // Seleccionando los 30 mejores atributos

    for (int i=0; i<30 and !mi_cola.empty(); i++){
      candidato = mi_cola.top();

      cout << candidato.nombre << ": " << candidato.valor << endl;
      atrib = candidato.atributo;

      variable_t nueva_variable(candidato.nombre,true,true,atrib);
      F.Add_Variable_Lista(nueva_variable,candidato.operacion,candidato.operando1,candidato.operando2,nueva_variable.Inf_Range(),nueva_variable.Sup_Range());

      mi_cola.pop();
    }

    //cin.get();

    delete [] py;
    delete [] Information;
    delete [] orden;
}




void Redefiniendo_Particiones_Variables (VectorVar &V, example_set &E)
{

    const char *intlabel[10]= {"2","3","4","5","6","7","8","9","10","11"};
    const char *intsensibilidad[4]= {"Uni","Izq","Cnt","Dch"};
    int n=V.TotalVariables(), conse=V.Consecuente();
    int ne;
    bool anadir;
    int ciclo;

    double *px;
    double *pz;
    double *py;
    double *pf;
    double **pxy;
    double **pzy;
    double **pfy;
    double Inf, H;
    double *Information, **entropia;
    int *orden;
    int **orden_clase;

    double measure[100];

    double minimo[V.SizeDomain(conse)],/* entropia[V.SizeDomain(conse)],*/ max_clase[V.SizeDomain(conse)], nuevo_max[V.SizeDomain(conse)];
    string NomVar1, NomVar2, NomVar3;
    domain_t atrib;

    int nv=V.N_Antecedente();
    Information = new double[nv];

    py = new double [V.SizeDomain(conse)];

    CalcularProbabilidadSimpleVariable(V,E,conse,V.SizeDomain(conse),py);

    double max=0, min=10000;
    for (int x=0; x<nv; x++)
    {
        px = new double [V.SizeDomain(x)];
        CalcularProbabilidadSimpleVariable(V,E,x,V.SizeDomain(x),px);
        pxy = new double*[V.SizeDomain(x)];
        for (int k=0; k<V.SizeDomain(x); k++)
            pxy[k]=new double [V.SizeDomain(conse)];
        CalcularProbabilidadConjunta2Variables(V,E,x,conse,V.SizeDomain(x),V.SizeDomain(conse),pxy);

        Information[x] = CalcularInformacion(V.SizeDomain(x), V.SizeDomain(conse), px, py, pxy, measure);
        //CalcularEntropia(V.SizeDomain(x), V.SizeDomain(conse), px, py, pxy, entropia);

        for (int k=0; k<V.SizeDomain(x); k++)
            delete [] pxy[k];
        delete [] pxy;
        delete [] px;

        //for (int i=0; i<V.SizeDomain(conse); i++)
        if (Information[x]>max)
            max=Information[x]; //+0.0001;

        NomVar1 = V.SPrintVar(x);

        cout << "Informacion (" << NomVar1 << ")= " << Information[x]<<"     ";
        for (int i=0; i<V.SizeDomain(conse); i++)
          cout << measure[i] <<" | ";

        cout << endl;
    }

    /*cout << "--------------------------------------------------------------------------------------------" << endl;
    for (int i=0; i<V.SizeDomain(conse); i++)
        cout << max_clase[i] << " | ";
    cout << endl;
    cout << "--------------------------------------------------------------------------------------------" << endl;*/

    // Trabajar solo con las 100 variables mas informadas
    /*cout << "--------------------------------------------------------------------------------------------" << endl;
    for (int i=0; i<nv && i<100; i++)
      cout << Information[orden[i]] << endl;
    cout << endl;
    cout << "--------------------------------------------------------------------------------------------" << endl;*/


    cout << "-------------------------------------  Buscando Particiones Mas Informadas -------------------------------------------------------" << endl;
    // Buscar particiones mas informadas
    int mejor_nlabels, mejor_sensi;
    double mejor_Inf;
    for (int x=0; x<nv; x++)
    {
        if(!V.IsDiscrete(x))
        {
            NomVar1= V.SPrintVar(x);
            mejor_nlabels=-1;
            mejor_Inf=Information[x];
            int op=15;
            for (int nlabels=5; nlabels<=5; nlabels++)
            {
                for (int sensibilidad=0; sensibilidad<4; sensibilidad++)
                {
                    atrib = Crear_Dominio(V,E,x,x,op,nlabels,sensibilidad);
                    pf = new double [atrib.Size()];
                    CalcularProbabilidadSimpleVariable_op(V,E,x,x,atrib,op,pf);
                    pfy = new double*[atrib.Size()];
                    for (int k=0; k<atrib.Size(); k++)
                        pfy[k]=new double [V.SizeDomain(conse)];
                    CalcularProbabilidadConjunta2Variables_op(V,E,x,x,atrib,conse,V.SizeDomain(conse),op, pfy);

                    Inf = CalcularInformacion(atrib.Size(), V.SizeDomain(conse), pf, py, pfy, measure);

                    if (mejor_Inf<Inf)
                    {
                        mejor_Inf=Inf;
                        mejor_nlabels=nlabels;
                        mejor_sensi=sensibilidad;
                    }
                    for (int k=0; k<atrib.Size(); k++)
                        delete [] pfy[k];
                    delete [] pfy;
                    delete [] pf;
                }
            }

            if (mejor_nlabels!=-1)
            {
                int unidad;
                double factor;

                unidad = V.Get_Unit(x);
                factor = V.Get_Convert_Factor(x);
                atrib = Crear_Dominio(V,E,x,x,op,mejor_nlabels,mejor_sensi);
                NomVar3 = NomVar1 + "_<"+intlabel[mejor_nlabels-2]+","+intsensibilidad[mejor_sensi]+">";
                cout << NomVar3 << endl;
                variable_t nueva_variable(NomVar1,true,true,atrib);
                V.Asigna(x,nueva_variable);
                V.Put_Unit(x,unidad);
                V.Put_Convert_Factor(x,factor);
                Information[x]=mejor_Inf;
            }
        }
    }
    delete [] Information;
    delete [] py;
}





void MedidasDeInformacion_operandoresUnarios(VectorVar &V, example_set &E, int nuevo,bool *bina)
{
    struct nova_var
    {
        int v1;
        int v2;
        int op;
    };

    const char *intlabel[10]= {"2","3","4","5","6","7","8","9","10","11"};
    const char *intsensibilidad[4]= {"Uni","Izq","Cnt","Dch"};
    int n=V.TotalVariables(), conse=V.Consecuente();
    int ne;
    bool anadir;
    int ciclo;

    double *px;
    double *pz;
    double *py;
    double *pf;
    double **pxy;
    double **pzy;
    double **pfy;
    double Inf, H;
    double *Information, **entropia;
    int *orden;
    int **orden_clase;

    double measure[100];

    double minimo[V.SizeDomain(conse)],/* entropia[V.SizeDomain(conse)],*/ max_clase[V.SizeDomain(conse)], nuevo_max[V.SizeDomain(conse)];
    nova_var a_introducir[V.SizeDomain(conse)];
    string NomVar1, NomVar2, NomVar3;
    domain_t atrib;

    for (int i=0; i<V.SizeDomain(conse); i++)
    {
        max_clase[i]=0;
        minimo[i]=100000.0;
        a_introducir[i].v1=-1;
        a_introducir[i].v2=-1;
        a_introducir[i].op=-1;
    }

    int nv=V.N_Antecedente();
    Information = new double[nv];
    entropia = new double*[V.SizeDomain(conse)];
    orden = new int[nv];
    orden_clase = new int*[V.SizeDomain(conse)];
    for (int i=0; i<V.SizeDomain(conse); i++)
    {
        orden_clase[i]= new int[nv];
        entropia[i]= new double [nv];
    }
    int nmax=0;

    py = new double [V.SizeDomain(conse)];

    CalcularProbabilidadSimpleVariable(V,E,conse,V.SizeDomain(conse),py);

    double max=0, min=10000;
    for (int x=0; x<nv; x++)
    {
        px = new double [V.SizeDomain(x)];
        CalcularProbabilidadSimpleVariable(V,E,x,V.SizeDomain(x),px);
        pxy = new double*[V.SizeDomain(x)];
        for (int k=0; k<V.SizeDomain(x); k++)
            pxy[k]=new double [V.SizeDomain(conse)];
        CalcularProbabilidadConjunta2Variables(V,E,x,conse,V.SizeDomain(x),V.SizeDomain(conse),pxy);

        Information[x] = CalcularInformacion(V.SizeDomain(x), V.SizeDomain(conse), px, py, pxy, measure);
        InsertarOrdenado (Information, orden, x);
        //CalcularEntropia(V.SizeDomain(x), V.SizeDomain(conse), px, py, pxy, entropia);

        for (int k=0; k<V.SizeDomain(x); k++)
            delete [] pxy[k];
        delete [] pxy;
        delete [] px;

        //for (int i=0; i<V.SizeDomain(conse); i++)
        if (Information[x]>max)
            max=Information[x]; //+0.0001;

        NomVar1 = V.SPrintVar(x);

        for (int i=0; i<V.SizeDomain(conse); i++)
        {
            entropia[i][x]=measure[i];
            InsertarOrdenado (entropia[i], orden_clase[i],x);
            //cout << entropia[i][x] << " | ";
            if (entropia[i][x]>max_clase[i])
            {
                max_clase[i]=entropia[i][x];
                nuevo_max[i]=max_clase[i];
            }
        }
        //cout << endl;
    }


    for (int x=0; x<nv; x++)
    {
        NomVar1 = V.SPrintVar(orden[x]);
        cout << "Informacion (" << NomVar1 << ")= " << Information[orden[x]]<<"     ";
        for (int i=0; i<V.SizeDomain(conse); i++)
        {
            cout << entropia[i][orden[x]] << " | ";
        }
        cout << endl;
    }
    cout << "--------------------------------------------------------------------------------------------" << endl;
    for (int i=0; i<V.SizeDomain(conse); i++)
        cout << max_clase[i] << " | ";
    cout << endl;
    cout << "--------------------------------------------------------------------------------------------" << endl;


    double mejor_max_clase[100];
    double mejor_Inf;
    int mejor_nlabels, mejor_sensi;
    cout << "-------------------------------------  Buscando Operadores Unarios -------------------------------------------------------" << endl;
    for (int z=nuevo; z<nv-1; z++)
    {
        //int z=orden[t1];
        if (!bina[z] && !V.IsDiscrete(z))
        {
            NomVar1= V.SPrintVar(z);
            cout << "-------------------------------------  " << NomVar1 << " -------------------------------------------------------" << endl;
            for (int op=6; op<12 /*9 */; op++)
            {
                mejor_nlabels=-1;
                mejor_Inf=Information[z];
                for (int i=0; i<V.SizeDomain(conse); i++)
                    mejor_max_clase[i]=max_clase[i];

                for (int nlabels=3; nlabels<4; nlabels++)
                {
                    for (int sensibilidad=0; sensibilidad<4; sensibilidad++)
                    {
                        atrib = Crear_Dominio(V,E,z,z,op,nlabels, sensibilidad);
                        pf = new double [atrib.Size()];
                        CalcularProbabilidadSimpleVariable_op(V,E,z,z,atrib,op,pf);
                        pfy = new double*[atrib.Size()];
                        for (int k=0; k<atrib.Size(); k++)
                            pfy[k]=new double [V.SizeDomain(conse)];
                        CalcularProbabilidadConjunta2Variables_op(V,E,z,z,atrib,conse,V.SizeDomain(conse),op, pfy);

                        Inf = CalcularInformacion(atrib.Size(), V.SizeDomain(conse), pf, py, pfy, measure);

                        for (int i=0; i<V.SizeDomain(conse) /*&& !mejora*/; i++)
                        {
                            if (measure[i]>mejor_max_clase[i]+0.001  && measure[i]<=1.0)
                            {
                                mejor_nlabels=nlabels;
                                mejor_sensi=sensibilidad;
                                mejor_max_clase[i]=measure[i];
                                cout << "mejora para la clase " << i << " (de " << max_clase[i] << " a " << measure[i] <<")\n";
                            }
                        }
                        for (int k=0; k<atrib.Size(); k++)
                            delete [] pfy[k];
                        delete [] pfy;
                        delete [] pf;
                    }
                }

                if (mejor_nlabels!=-1  && conse<999)
                {
                    atrib = Crear_Dominio(V,E,z,z,op,mejor_nlabels,mejor_sensi);
                    pf = new double [atrib.Size()];
                    CalcularProbabilidadSimpleVariable_op(V,E,z,z,atrib,op,pf);
                    pfy = new double*[atrib.Size()];
                    for (int k=0; k<atrib.Size(); k++)
                        pfy[k]=new double [V.SizeDomain(conse)];
                    CalcularProbabilidadConjunta2Variables_op(V,E,z,z,atrib,conse,V.SizeDomain(conse),op, pfy);

                    Inf = CalcularInformacion(atrib.Size(), V.SizeDomain(conse), pf, py, pfy, measure);
                    nmax++;
                    NomVar1 = V.SPrintVar(z);
                    NomVar3 = Simbolo_Operador(op)+"("+NomVar1+")";
                    cout << NomVar3 << "<" << mejor_nlabels <<","<< intsensibilidad[mejor_sensi] <<">" << endl;
                    variable_t nueva_variable(NomVar3,true,true,atrib);
                    //nueva_variable.Pinta();
                    cout << endl;
                    if (op==6 || op==13)
                    {
                        double min,max;
                        Max_Min(E,V,z,z,-1,min,max);
                        V.Add_Variable(nueva_variable,op,z,z,min,max);
                    }
                    else
                        V.Add_Variable(nueva_variable,op,z,z,0,0);
                    E.AddNewColumn(z,z,op);
                    bina[conse]=false;
                    conse++; //despues de meter la nueva variable, conse se incrementa en una columna
                    for (int i=0; i<V.SizeDomain(conse); i++)
                    {
                        if (measure[i]>max_clase[i])
                        {
                            nuevo_max[i]=measure[i];
                            max_clase[i]=measure[i];  // para podar la busqueda
                            mejor_max_clase[i]=max_clase[i];
                        }
                    }
                    for (int k=0; k<atrib.Size(); k++)
                        delete [] pfy[k];
                    delete [] pfy;
                    delete [] pf;
                }
            }
        }
    }

    cout << "--------------------------------------------------------------------------------------------" << endl;
    for (int i=0; i<V.SizeDomain(conse); i++)
        cout << max_clase[i] << " | ";
    cout << endl;
    cout << "--------------------------------------------------------------------------------------------" << endl;

    for (int i=0; i<V.SizeDomain(conse); i++)
        delete [] orden_clase[i];
    delete [] orden_clase;

    delete [] py;
    delete [] Information;
    delete [] orden;

}





// Calculo de la ganancia de informacion de nuevos atributos creados a partir de operadores
void MedidasDeInformacion_operandores2(VectorVar &V, example_set &E, int nuevo,bool *bina)
{
    struct nova_var
    {
        int v1;
        int v2;
        int op;
    };

    const char *intlabel[10]= {"2","3","4","5","6","7","8","9","10","11"};
    const char *intsensibilidad[4]= {"Uni","Izq","Cnt","Dch"};
    int n=V.TotalVariables(), conse=V.Consecuente();
    int ne;
    bool anadir;
    int ciclo;

    double *px;
    double *pz;
    double *py;
    double *pf;
    double **pxy;
    double **pzy;
    double **pfy;
    double Inf, H;
    double *Information, **entropia;
    int *orden;
    int **orden_clase;

    double measure[100];

    double minimo[V.SizeDomain(conse)],/* entropia[V.SizeDomain(conse)],*/ max_clase[V.SizeDomain(conse)], nuevo_max[V.SizeDomain(conse)];
    nova_var a_introducir[V.SizeDomain(conse)];
    string NomVar1, NomVar2, NomVar3;
    domain_t atrib;

    for (int i=0; i<V.SizeDomain(conse); i++)
    {
        max_clase[i]=0;
        minimo[i]=100000.0;
        a_introducir[i].v1=-1;
        a_introducir[i].v2=-1;
        a_introducir[i].op=-1;
    }

    int nv=V.N_Antecedente();
    Information = new double[nv];
    entropia = new double*[V.SizeDomain(conse)];
    orden = new int[nv];
    orden_clase = new int*[V.SizeDomain(conse)];
    for (int i=0; i<V.SizeDomain(conse); i++)
    {
        orden_clase[i]= new int[nv];
        entropia[i]= new double [nv];
    }
    int nmax=0;

    py = new double [V.SizeDomain(conse)];

    CalcularProbabilidadSimpleVariable(V,E,conse,V.SizeDomain(conse),py);

    double max=0, min=10000;
    for (int x=0; x<nv; x++)
    {
        px = new double [V.SizeDomain(x)];
        CalcularProbabilidadSimpleVariable(V,E,x,V.SizeDomain(x),px);
        pxy = new double*[V.SizeDomain(x)];
        for (int k=0; k<V.SizeDomain(x); k++)
            pxy[k]=new double [V.SizeDomain(conse)];
        CalcularProbabilidadConjunta2Variables(V,E,x,conse,V.SizeDomain(x),V.SizeDomain(conse),pxy);

        Information[x] = CalcularInformacion(V.SizeDomain(x), V.SizeDomain(conse), px, py, pxy, measure);
        InsertarOrdenado (Information, orden, x);
        //CalcularEntropia(V.SizeDomain(x), V.SizeDomain(conse), px, py, pxy, entropia);

        for (int k=0; k<V.SizeDomain(x); k++)
            delete [] pxy[k];
        delete [] pxy;
        delete [] px;

        //for (int i=0; i<V.SizeDomain(conse); i++)
        if (Information[x]>max)
            max=Information[x]; //+0.0001;

        NomVar1 = V.SPrintVar(x);

        for (int i=0; i<V.SizeDomain(conse); i++)
        {
            entropia[i][x]=measure[i];
            InsertarOrdenado (entropia[i], orden_clase[i],x);
            //cout << entropia[i][x] << " | ";
            if (entropia[i][x]>max_clase[i])
            {
                max_clase[i]=entropia[i][x];
                nuevo_max[i]=max_clase[i];
            }
        }
        //cout << endl;
    }


    for (int x=0; x<nv; x++)
    {
        NomVar1 = V.SPrintVar(orden[x]);
        cout << "Informacion (" << NomVar1 << ")= " << Information[orden[x]]<<"     ";
        for (int i=0; i<V.SizeDomain(conse); i++)
        {
            cout << entropia[i][orden[x]] << " | ";
        }
        cout << endl;
    }
    cout << "--------------------------------------------------------------------------------------------" << endl;
    for (int i=0; i<V.SizeDomain(conse); i++)
    {
        cout << max_clase[i] << " | ";
        nuevo_max[i]=max_clase[i];
    }
    cout << endl;
    cout << "--------------------------------------------------------------------------------------------" << endl;


    double mejor_max_clase[100];
    double mejor_Inf;
    int mejor_nlabels, mejor_sensi;
    cout << "-------------------------------------  Buscando Condiciones -------------------------------------------------------" << endl;
    for (int x=0; x<nv-1; x++)
    {
        ///int x=orden[t1];
        NomVar1= V.SPrintVar(x);
        if (!bina[x] && !V.IsDiscrete(x))
        {
            cout <<"----------------------------------------------------  " << NomVar1 <<" ----------------------------\n";

            for (int z=nuevo; z<nv; z++)
            {
                //int z=orden[z];
                if (!bina[z] && !V.IsDiscrete(z))
                {
                    for (int op=1; op<2; op++)
                    {
                        atrib = Crear_Dominio(V,E,x,z,op);
                        pf = new double [atrib.Size()];
                        CalcularProbabilidadSimpleVariable_op(V,E,x,z,atrib,op,pf);
                        pfy = new double*[atrib.Size()];
                        for (int k=0; k<atrib.Size(); k++)
                            pfy[k]=new double [V.SizeDomain(conse)];
                        CalcularProbabilidadConjunta2Variables_op(V,E,x,z,atrib,conse,V.SizeDomain(conse),op, pfy);

                        Inf = CalcularInformacion(atrib.Size(), V.SizeDomain(conse), pf, py, pfy, measure);
                        //CalcularEntropia(atrib.Size(), V.SizeDomain(conse), pf, py, pfy, entropia);

                        bool mejora=false;
                        for (int i=0; i<V.SizeDomain(conse) /*&& !mejora*/; i++)
                        {
                            if (measure[i]>max_clase[i]+0.001+99*(nuevo_max[i]-max_clase[i])/100.0)
                            {
                                //max_clase[i]=entropia[i];
                                mejora = true;
                                cout << "mejora para la clase " << i << " (de " << max_clase[i]+0.001+99*(nuevo_max[i]-max_clase[i])/100.0 << " a " << measure[i] <<")\n";
                                //cout << "Nuevo valor: " << entropia[i] << "   Anterior Valor: " << max_clase[i] << endl;
                            }
                        }

                        if (mejora && conse<999)
                        {
                            nmax++;
                            NomVar1 = V.SPrintVar(x);
                            NomVar2 = V.SPrintVar(z);
                            NomVar3 = "("+NomVar1+Simbolo_Operador(op)+NomVar2+")";
                            cout << NomVar3 << endl;
                            variable_t nueva_variable(NomVar3,true,true,atrib);
                            V.Add_Variable(nueva_variable,op,x,z,0,1);
                            E.AddNewColumn(x,z,op);
                            if (op<=1)
                                bina[conse]=true;
                            else
                                bina[conse]=bina[x] || bina[z];
                            conse++; //despues de meter la nueva variable, conse se incrementa en una columna
                            //max=Inf+0.0001;
                            //max=Inf;
                            for (int i=0; i<V.SizeDomain(conse); i++)
                                if (measure[i]>max_clase[i] /*&& Inf<=1.0*/)
                                {
                                    if (measure[i]>nuevo_max[i])
                                        nuevo_max[i]=measure[i];
                                    //max_clase[i]=measure[i];  // para podar la busqueda
                                }
                        }


                        for (int k=0; k<atrib.Size(); k++)
                            delete [] pfy[k];
                        delete [] pfy;
                        delete [] pf;

                    }
                }
            }
        }
    }
    cout << "-------------------------------------  Uniendo Condiciones con Operadores Logicos -------------------------------------------------------" << endl;
    for (int i=0; i<V.SizeDomain(conse); i++)
    {
        cout << max_clase[i] << " | ";
        max_clase[i]=nuevo_max[i];
    }
    cout << endl;

    for (int x=0; x<nv-1; x++)
    {
        ///int x=orden[t1];
        if (bina[x])
        {
            NomVar1= V.SPrintVar(x);
            cout <<"----------------------------------------------------  " << NomVar1 <<" ----------------------------\n";

            if (x>=nuevo)
                nuevo=x+1;
            for (int z=nuevo; z<nv; z++)
            {
                //int z=orden[z];
                if (bina[z])
                {
                    for (int op=20; op<22; op++)
                    {
                        atrib = Crear_Dominio(V,E,x,z,op);
                        pf = new double [atrib.Size()];
                        CalcularProbabilidadSimpleVariable_op(V,E,x,z,atrib,op,pf);
                        pfy = new double*[atrib.Size()];
                        for (int k=0; k<atrib.Size(); k++)
                            pfy[k]=new double [V.SizeDomain(conse)];
                        CalcularProbabilidadConjunta2Variables_op(V,E,x,z,atrib,conse,V.SizeDomain(conse),op, pfy);

                        Inf = CalcularInformacion(atrib.Size(), V.SizeDomain(conse), pf, py, pfy, measure);
                        //CalcularEntropia(atrib.Size(), V.SizeDomain(conse), pf, py, pfy, entropia);

                        bool mejora=false;
                        for (int i=0; i<V.SizeDomain(conse) /*&& !mejora*/; i++)
                        {
                            if (measure[i]>max_clase[i]+0.001+99*(nuevo_max[i]-max_clase[i])/100.0)
                            {
                                //max_clase[i]=entropia[i];
                                mejora = true;
                                cout << "mejora para la clase " << i << " (de " << max_clase[i]+0.001+99*(nuevo_max[i]-max_clase[i])/100.0 << " a " << measure[i] <<")\n";
                                //cout << "Nuevo valor: " << entropia[i] << "   Anterior Valor: " << max_clase[i] << endl;
                            }
                        }

                        if (mejora && conse<999)
                        {
                            nmax++;
                            NomVar1 = V.SPrintVar(x);
                            NomVar2 = V.SPrintVar(z);
                            NomVar3 = "("+NomVar1+Simbolo_Operador(op)+NomVar2+")";
                            cout << NomVar3 << endl;
                            variable_t nueva_variable(NomVar3,true,true,atrib);
                            V.Add_Variable(nueva_variable,op,x,z,0,0);
                            E.AddNewColumn(x,z,op);
                            if (op<=1)
                                bina[conse]=true;
                            else
                                bina[conse]=bina[x] || bina[z];
                            conse++; //despues de meter la nueva variable, conse se incrementa en una columna
                            //max=Inf+0.0001;
                            //max=Inf;
                            for (int i=0; i<V.SizeDomain(conse); i++)
                                if (measure[i]>max_clase[i] /*&& Inf<=1.0*/)
                                {
                                    nuevo_max[i]=measure[i];
                                    //max_clase[i]=measure[i];  // para podar la busqueda
                                }
                        }


                        for (int k=0; k<atrib.Size(); k++)
                            delete [] pfy[k];
                        delete [] pfy;
                        delete [] pf;

                    }
                }
            }
        }
    }



    /* cout << "-------------------------------------  Buscando Funciones Binarias -------------------------------------------------------" << endl;
    int fin,x;
    if (nuevo==0){
        fin=nuevo+20; //5 para ImagenJuanCa
        if (fin>nv-1)
      fin=nv-1;
    }
    else
        fin=nv-1;

    for (int t=fin; t>=nuevo; t--){
      if (nuevo==0)
         x=orden[t];
      else
         x=t;
      NomVar1= V.SPrintVar(x);
      if (!bina[x] && !V.IsDiscrete(x)){
        cout <<"----------------------------------------------------  " << NomVar1 << " ("<<Information[x] <<") ----------------------------\n";

       for (int z=fin; z<nv; z++){
          //int z=orden[t2];
          if (!bina[z]&& !V.IsDiscrete(z) && z!=x){
            for (int i=0; i<V.SizeDomain(conse); i++)
         mejor_max_clase[i]=max_clase[i];

            for (int op=2; op<4; op++){
        mejor_nlabels=-1;
        mejor_Inf=Information[z];
              for (int nlabels=3; nlabels<4; nlabels++){
           for (int sensibilidad=0; sensibilidad<1; sensibilidad++){
                      atrib = Crear_Dominio(V,E,x,z,op,nlabels,sensibilidad);
                      pf = new double [atrib.Size()];
                     CalcularProbabilidadSimpleVariable_op(V,E,x,z,atrib,op,pf);
                     pfy = new double*[atrib.Size()];
                     for (int k=0; k<atrib.Size(); k++)
                         pfy[k]=new double [V.SizeDomain(conse)];
                     CalcularProbabilidadConjunta2Variables_op(V,E,x,z,atrib,conse,V.SizeDomain(conse),op, pfy);

                     Inf = CalcularInformacion(atrib.Size(), V.SizeDomain(conse), pf, py, pfy, measure);


              for (int i=0; i<V.SizeDomain(conse); i++){
                       if (measure[i]>mejor_max_clase[i]+0.001){
    	       mejor_nlabels=nlabels;
    	       mejor_sensi = sensibilidad;
    	       mejor_max_clase[i]=measure[i];
                   cout << "mejora para la clase " << i << " (de " << max_clase[i] << " a " << measure[i] <<")\n";
                }
                    }
               for (int k=0; k<atrib.Size(); k++)
                  delete [] pfy[k];
              delete [] pfy;
              delete [] pf;
        }
    }

          if (mejor_nlabels!=-1  && conse<999){
             atrib = Crear_Dominio(V,E,x,z,op,mejor_nlabels,mejor_sensi);
             pf = new double [atrib.Size()];
             CalcularProbabilidadSimpleVariable_op(V,E,x,z,atrib,op,pf);
             pfy = new double*[atrib.Size()];
             for (int k=0; k<atrib.Size(); k++)
                  pfy[k]=new double [V.SizeDomain(conse)];
             CalcularProbabilidadConjunta2Variables_op(V,E,x,z,atrib,conse,V.SizeDomain(conse),op, pfy);

             Inf = CalcularInformacion(atrib.Size(), V.SizeDomain(conse), pf, py, pfy, measure);
       nmax++;
             NomVar1 = V.SPrintVar(x);
             NomVar2 = V.SPrintVar(z);
             NomVar3 = "("+NomVar1+Simbolo_Operador(op)+NomVar2+")";
        cout << NomVar3 << endl;
                variable_t nueva_variable(NomVar3,true,true,atrib);
                V.Add_Variable(nueva_variable,op,x,z,0,1);
                E.AddNewColumn(x,z,op);
          if (op<=1)
            bina[conse]=true;
         else
           bina[conse]=bina[x] || bina[z];
               conse++; //despues de meter la nueva variable, conse se incrementa en una columna

              for (int i=0; i<V.SizeDomain(conse); i++){
          if (measure[i]>max_clase[i]){
            nuevo_max[i]=measure[i];
            max_clase[i]=measure[i];  // para podar la busqueda
            mejor_max_clase[i]=max_clase[i];
         }
        }
             for (int k=0; k<atrib.Size(); k++)
               delete [] pfy[k];
             delete [] pfy;
             delete [] pf;
          }
        }
       }
      }
     }
    } */





    cout << "Numero de atributos prometedores= " << nmax << endl;
    for (int i=0; i<V.SizeDomain(conse); i++)
        delete [] orden_clase[i];
    delete [] orden_clase;

    delete [] py;
    delete [] Information;
    delete [] orden;
}




void MedidasDeInformacion_2variables_cond_consecuente(VectorVar &V, example_set &E)
{

    int conse=V.Consecuente();
    int nv=V.N_Antecedente();
    /*Calculo la distribucion de probabilidad de la clase */
    double *Pz = new double [V.SizeDomain(conse)], *Px, *Py;
    double **Pxz, **Pyz;
    double ***Pxyz;

    double Info, *Info_z=new double[V.SizeDomain(conse)], Info_x, Info_y;

    string NomX, NomY, NomZ;

    NomZ = V.SPrintVar(conse);

    CalcularProbabilidadSimpleVariable(V,E,conse,V.SizeDomain(conse),Pz);

    /*Calculo las probabilidades conjuntas condicionadas a la clase de todas las variables */
    for (int x=0; x<nv; x++)
    {
        NomX = V.SPrintVar(x);
        Px = new double[V.SizeDomain(x)];
        CalcularProbabilidadSimpleVariable(V,E,x,V.SizeDomain(x),Px);

        Pxz = new double*[V.SizeDomain(x)];
        for (int i=0; i<V.SizeDomain(x); i++)
            Pxz[i]=new double[V.SizeDomain(conse)];

        CalcularProbabilidadConjunta2Variables(V,E,x,conse,V.SizeDomain(x),V.SizeDomain(conse),Pxz);
        Info_x = CalcularInformacion (V.SizeDomain(x), V.SizeDomain(conse), Px, Pz, Pxz);

        for (int y=x+1; y<nv; y++)
        {
            NomY = V.SPrintVar(y);
            Py = new double [V.SizeDomain(y)];
            CalcularProbabilidadSimpleVariable(V,E,y,V.SizeDomain(y),Py);

            Pyz = new double*[V.SizeDomain(y)];
            for (int i=0; i<V.SizeDomain(y); i++)
                Pyz[i]=new double[V.SizeDomain(conse)];

            CalcularProbabilidadConjunta2Variables(V,E,y,conse,V.SizeDomain(y),V.SizeDomain(conse),Pyz);
            Info_y = CalcularInformacion (V.SizeDomain(y), V.SizeDomain(conse), Py, Pz, Pyz);

            Pxyz= new double **[V.SizeDomain(x)];
            for (int i=0; i<V.SizeDomain(x); i++)
            {
                Pxyz[i]=new double*[V.SizeDomain(y)];
                for (int j=0; j<V.SizeDomain(y); j++)
                    Pxyz[i][j]= new double[V.SizeDomain(conse)];
            }

            CalcularProbabilidadConjunta3Variables(V,E,x,y,conse, V.SizeDomain(x), V.SizeDomain(y), V.SizeDomain(conse), Pxyz);

            Info = CalcularInformacionCondicionada2variablesRespectoClase (V, E, V.SizeDomain(x), V.SizeDomain(y), V.SizeDomain(conse), Pxyz, Pxz, Pyz,Pz,Info_z);

            if (Info >= Info_x && Info >=Info_y)
            {
                cout << "I("<<NomX<<";"<<NomY<<"/"<<NomZ<<")= " << Info << "    [";
                cout << Info_z[0];
                for (int i=1; i<V.SizeDomain(conse); i++)
                    cout << "|" <<Info_z[i];

                cout << "] " << "Info("<<NomX<<")= " << Info_x << "         Info("<<NomY<<")= " << Info_y << "    Suma= " << Info_x+Info_y << endl;
            }
            /* Eliminar memoria de Pyz */
            for (int i=0; i<V.SizeDomain(y); i++)
                delete [] Pyz[i];
            delete [] Pyz;

            /* Eliminar memoria de Pxyz */
            for (int i=0; i<V.SizeDomain(x); i++)
            {
                for (int j=0; j<V.SizeDomain(y); j++)
                    delete [] Pxyz[i][j];
                delete Pxyz[i];
            }
            delete [] Pxyz;
            delete [] Py;

        }
        for (int k=0; k<V.SizeDomain(x); k++)
            delete [] Pxz[k];
        delete [] Pxz;
        delete [] Px;
    }
    delete [] Pz;
    delete [] Info_z;
}





double ***Crear_Tabla_Adaptaciones(const VectorVar &V, const example_set &E)
{
    int n=V.TotalVariables();
    int m=E.N_Examples();
    int l;
    double valor;

    // Reservando la Memoria para las adaptaciones
    double ***Tabla = new double**[n];
    for (int i=0; i<n; i++)
    {
        l=V.SizeDomain(i);
        Tabla[i]= new double *[l];
        for (int j=0; j<l; j++)
            Tabla[i][j]= new double [m];
    }


    // Calculando los valores de adaptaciones de cada ejemplo a cada etiqueta de cada variable
    int valores_distinto_cero;
    for (int e=0; e<m; e++)
    {
        for (int v=0; v<n; v++)
        {
            valor=E.Data(e,v);
            l=V.SizeDomain(v);
            valores_distinto_cero=0;
            for (int j=0; j<l; j++){
              Tabla[v][j][e]=V.Adaptacion(valor,v,j);
              if (Tabla[v][j][e]>0 and Tabla[v][j][e]<=1)
                valores_distinto_cero++;
            }
            if (valores_distinto_cero==0){
              cout << "ERROR!!!\n";
              V.PrintVar(v);
              cout << "\nEl ejemplo " << e << " de la variable " << v << " de valor " << valor << " no tiene ninguna adaptacion con su dominio\n";
              exit(1);
            }

        }
    }

    return Tabla;
}


double ***Crear_Tabla_Adaptaciones_Nuevas(const VectorVar &V, const example_set &E)
{
    int n=V.TotalVariables();
    int m=E.N_Examples();
    int l;
    int var1, var2, op;
    double valor1, valor2, resultado;

    // Reservando la Memoria para las adaptaciones
    double ***Tabla = new double**[n];
    for (int i=0; i<n; i++)
    {
        l=V.SizeDomain(i);
        //V.Variable(i).Pinta();
        Tabla[i]= new double *[l];
        for (int j=0; j<l; j++){
            Tabla[i][j]= new double [m];
        }
    }


    // Calculando los valores de adaptaciones de cada ejemplo a cada etiqueta de cada variable
    int valores_distinto_cero;
    for (int e=0; e<m; e++)
    {
        for (int v=0; v<n; v++)
        {
            V.Funcion(v,var1,var2,op);
            l=V.SizeDomain(v);
            if (op==5 and E.Data(e,var2)==0)
              for (int j=0; j<l; j++)
                 Tabla[v][j][e]=0;
            else {
              valores_distinto_cero=0;
              resultado=Operacion(E.Data(e,var1),E.Data(e,var2), op);
              for (int j=0; j<l; j++){
                Tabla[v][j][e]=V.Adaptacion(resultado,v,j);
                if (Tabla[v][j][e]>0 and Tabla[v][j][e]<=1)
                  valores_distinto_cero++;
              }
              if (valores_distinto_cero==0){
                cout << "ERROR!!!\n";
                V.PrintVar(v);
                cout << "\nEl ejemplo " << e << " de la variable " << v << " de valor " << resultado << " no tiene ninguna adaptacion\n";
              exit(1);
            }
            }
        }
    }

    return Tabla;
}


void Destruir_Tabla_Adaptaciones(const VectorVar &V, const example_set &E, double ***&Tabla)
{
    int n=V.TotalVariables();
    int l;
    for (int i=0; i<n; i++)
    {
        l=V.SizeDomain(i);
        for (int j=0; j<l; j++)
            delete [] Tabla[i][j];
        delete [] Tabla[i];
    }
    delete [] Tabla;
}

//Adaptacion de un ejemplo con una etiqueta de la variable
double T_Adaptacion_Variable_Etiqueta(const VectorVar &V, double ***Tabla, int ejemplo, int variable, int etiqueta)
{
    return Tabla[variable][etiqueta][ejemplo];
}




//Adaptacion de un ejemplo con un variable
double T_Adaptacion_Variable(const VectorVar &V, double ***Tabla, int ejemplo, int variable, string labels)
{
    int l=V.SizeDomain(variable);
    double max=0;
    for (int etiqueta=0; etiqueta<l && max<1; etiqueta++)
        if ( labels[etiqueta]=='1' )
          if (max==0)
            max = Tabla[variable][etiqueta][ejemplo];
          else if (Tabla[variable][etiqueta][ejemplo]>0)
                  max=1;
    return max;
}



//Adaptacion de un ejemplo al antecedente de una regla considerando todas las variables relevantes
double T_Adaptacion_Antecedente(const VectorVar &V, double ***Tabla, int ejemplo, string regla, double umbral)
{
    int n=V.N_Antecedente();
    double max=1,aux;
    string sub;
    int trozo=0,tam, unos;
    for (int v=0; v<n && max>0; v++)
    {
        tam=V.SizeDomain(v);
        sub=regla.substr(trozo,tam);
        aux = T_Adaptacion_Variable(V,Tabla,ejemplo,v,sub);
        if (aux<max)
            max=aux;
        trozo+=tam;
    }
    return max;
}




//Adaptacion de un ejemplo al antecedente de una regla
double T_Adaptacion_Antecedente(const VectorVar &V, double ***Tabla, int ejemplo, string regla,
                                double *nivel_variable, double umbral)
{
    int n=V.N_Antecedente();
    double max=1,aux;
    string sub;
    int trozo=0,tam, unos;
    for (int v=0; v<n && max>0; v++)
    {
        tam=V.SizeDomain(v);
        if (nivel_variable[v]>=umbral)
        {
            sub=regla.substr(trozo,tam);
            aux = T_Adaptacion_Variable(V,Tabla,ejemplo,v,sub);
            if (aux<max)
                max=aux;
        }
        trozo+=tam;
    }
    return max;
}


//Adaptacion de un ejemplo al antecedente de una regla
double T_Adaptacion_Antecedente_Nuevas(const VectorVar &V, double ***Tabla, int ejemplo, string regla,
                                       int *nivel_variable_funciones, int tama_nivel_variable_funciones, double matching_previo)
{
    int n=V.N_Antecedente();
    double max=matching_previo,aux;
    string sub;
    int trozo=0,tam, unos;
    for (int func=0; (max>0) && (func<tama_nivel_variable_funciones); func++)
    {
        if (nivel_variable_funciones[func]!=0)
        {
            tam=V.SizeDomain(nivel_variable_funciones[func]-1);
            sub=regla.substr(trozo,tam);
            aux = T_Adaptacion_Variable(V,Tabla,ejemplo,nivel_variable_funciones[func]-1,sub);
            if (aux<max)
                max=aux;
        }
        //else
           tam=2;

        trozo+=tam;
    }
    return max;
}



//Adaptacion de un ejemplo con el consecuente
void T_AdaptacionC(const VectorVar &V, double ***Tabla, int ejemplo, int etiq, double &pos, double &neg)
{
    int conse=V.Consecuente();
    pos =Tabla[conse][etiq][ejemplo];
    if (pos==1)
        neg=0;
    else
        neg=1;
}


//Mejor etiqueta dada una variable y un ejemplo
int BetterEtiqueta(const VectorVar &V,double ***Tabla, int variable, int ejemplo){

    int l=V.SizeDomain(variable), et=0;
    double max=0;

    for (int etiqueta=0; etiqueta<l && max<1; etiqueta++){
        if (Tabla[variable][etiqueta][ejemplo]>max){
            max = Tabla[variable][etiqueta][ejemplo];
            et = etiqueta;
        }
    }
    return et;
}









double Bondad_old (multipoblacion &G, const VectorVar &V, const example_set &E, ruleset &R, double *valor, int N_valor,
               int elemento_i, bool discriminante, bool reducir_reglas, int *ponderacion_por_clase, int nrule,
               double ***Tabla, double ***Tabla_Nuevas, VectorVar &F, int *basica_para_regla, int *errores_para_regla)
{

    int b1, b2, n1, n2, r1, r2, r3;
    char *nb1, *nb2;
    int *nn1,*nn2,r;
    double *nr1, *nr2, *nr3;
    string regla, regla_funciones;
    int var1, var2;
    int ne=E.N_Examples();
    vectordouble w(V.N_Antecedente());
    double positivos=0, negativos=0, rango, eficacia, negativos_que_conllevan_error=0, positivos_a_incrementar=0,
                                  positivos_sobrecubrimiento=0, negativos2=0;
    int j;
    double aciertos=0, fallos=0, aciertos_aqui=0, fallos_aqui=0, sin_cubrir=0, aciertos_no_cubiertos=0, max_pos=0,
                              max_neg=0, max_neg2=0;
    int clase = G.ClaseIndividuo(elemento_i);
    bool menos_infinito=true, valida=true;
    int variables_irrelevantes=0, funciones_irrelevantes=0, variables_estables=0, funciones_estables=0, ceros_variables=0, ceros_funciones=0, unos_variables=0, unos_funciones=0, variables_relevantes=0, funciones_relevantes=0;


    nb1=G.Subpoblacion_Binaria(0,elemento_i,b1);
    char *s= new char[b1+1];
    valor[2]=0;  //numero de unos en el cromosoma
    for (j=0; j<b1; j++)
    {
        s[j]=nb1[j];
        if (nb1[j]=='1')
            valor[2]++;
        s[b1]='\0';
    }
    regla=s;
    delete []s;


    nb2=G.Subpoblacion_Binaria(1,elemento_i,b2);
    s= new char[b2+1];
    for (j=0; j<b2; j++)
    {
        s[j]=nb2[j];
        if (nb2[j]=='1')
            valor[2]++;
        s[b2]='\0';
    }
    regla_funciones=s;
    delete []s;


    nn1=G.Subpoblacion_Entera(0,elemento_i,n1);
    nn2=G.Subpoblacion_Entera(1,elemento_i,n2);
    nr1=G.Subpoblacion_Real(0,elemento_i,r1);
    nr2 = new double[ne];


    valor[1]=0; //numero de variable irrelevantes
    valor[2]=0; //funciones irrelevantes
    valor[3]=0; //asignaciones comprensibles
    valor[4]=0; //numero de ceros
    int pos=0; // posicion en el cromosoma nivel de valor
    int tamanio; // numero de genes de una variable
    for (j=0; j<r1-1; j++)
    {
        tamanio= V.Variable(j).N_etiquetas();
        if (nr1[j]<nr1[r1-1])
            variables_irrelevantes++;
        else
        {
           if (Relevante(&(nb1[pos]),tamanio)){
              if (Asignacion_Comprensible(&(nb1[pos]),tamanio)){
                variables_estables++;
                ceros_variables+=NumeroCeros(&(nb1[pos]),tamanio);
                unos_variables+=NumeroUnos(&(nb1[pos]),tamanio);
              }
              else {
                // Asignacion no comprensible
                // considero entonces la regla como invalida
                // valida=false;
              }
           }
           else {
              nr1[j]=nr1[r1-1]-0.001;
              variables_irrelevantes++;
           }

        }
        pos=pos+tamanio;
    }

    variables_relevantes = (r1-1) - variables_irrelevantes;

    // Quitar funciones repetidas

    for (int i=1; i<n2; i++){
      if (nn2[i]!=0) {
        int t=0;
        while (t<i and nn2[t]!=nn2[i])
          t++;
        if (t<i)
          nn2[i]=0;
      }
    }


    pos=0;
    for (j=0; j<n2; j++)
    {
        if (nn2[j]!=0)
        {
           tamanio= F.SizeDomain(nn2[j]-1);
           if (Relevante(&(nb2[pos]),tamanio)){
              if (Asignacion_Comprensible(&(nb2[pos]),tamanio)){
                 funciones_estables++;
                 ceros_funciones+=NumeroCeros(&(nb2[pos]),tamanio);
                 unos_funciones+=NumeroUnos(&(nb2[pos]),tamanio);
              }
              else {
                // Asignacion no comprensible
                // considero entonces la regla como invalida
                // valida=false;
              }
           }
           else {
              funciones_irrelevantes++;
              nn2[j]=0;
           }
        }
        else {
            tamanio=2;
            funciones_irrelevantes++;
        }
        tamanio=2;
        pos=pos+tamanio;
    }

    funciones_relevantes = n2 - funciones_irrelevantes;


//Distingue entre descripciones discriminantes o descriptivas
    /*if (!discriminante)
        valor[1]=r1-1+n2-valor[1]-valor[2];*/



    int *basica_esta_regla = 0;
    int *errores_esta_regla = 0;
    if (nrule>0) {
      basica_esta_regla = new int[nrule];
      errores_esta_regla = new int[nrule];
      for (int i=0; i<nrule; i++){
        basica_esta_regla[i]=0;
        errores_esta_regla[i]=0;
      }
    }

    int *ejemplos_afectado_por_esta_regla = new int[ne];

//double *aux = new double[ne];
    double aux_p;
    double aux_n;
    double peso, kk;
    double fitness0_global_regla=-ne;
// Comprobar que es una regla valida
    bool esta_cubierto;
//bool regla_valida = V.Es_Valida(regla,nr1,nr1[r1-1],valor[2]);
    bool regla_valida = V.Es_Valida(regla,nr1,nr1[r1-1],kk);
    double afecta=0;

    if (valida)
       valida = (ponderacion_por_clase[clase]!=0);

// Calculo previo para obtener el peso de la regla sobre atributos originales
    if (regla_valida and valida)
    {
        for (j=0; j<ne; j++)
        {

            T_AdaptacionC(V,Tabla,j,nn1[0],aux_p,aux_n);

            //esta_cubierto = E.Is_Covered(j);

            nr2[j]=T_Adaptacion_Antecedente(V,Tabla,j,regla,nr1,nr1[r1-1]);
            //nr2[j]=1;
            nr2[j]=T_Adaptacion_Antecedente_Nuevas(F, Tabla_Nuevas, j, regla_funciones, nn2, n2, nr2[j]);

            if (nr2[j]*aux_p>0)
              afecta++;

            positivos=positivos+(nr2[j]*aux_p);
            negativos=negativos+(nr2[j]*aux_n);
            if (aux_p<aux_n)
                nr2[j]=-nr2[j];

        }

        //for (int i=0; i<nrule; i++)
          //cout <<"Regla " << i <<": " << basica_para_regla[i] << endl;
        //fitness0_global_regla = (positivos_global-negativos_global)/(ponderacion_por_clase[clase]);
    }


//peso = (positivos+1) / (positivos+negativos+1);
//peso = 1.0;
    peso = (positivos+1) / (positivos+negativos+1);
//peso = nr3[0];

 if (positivos>0)
   menos_infinito=false;


 /*if (afecta< 0.20*ponderacion_por_clase[clase])
   valida=false;*/
 if (ponderacion_por_clase[clase]>0)
   afecta = afecta / ponderacion_por_clase[clase];
 else
   valida=false;

 double positivos_global=0, negativos_global=0;
// Dado el peso de la regla obtener la bondad de la misma.
    int positivos_aux=0;
    if (regla_valida and valida)
    {
        positivos=0;
        negativos=0;
        for (j=0; j<ne; j++)
        {
            if (nr2[j]>=0)
              positivos_global+= (nr2[j]*peso);
            else
              negativos_global+= (-nr2[j]*peso);

            ejemplos_afectado_por_esta_regla[j]=0;

            esta_cubierto = E.Is_Covered(j);
            if (!esta_cubierto)
            {
                if (nr2[j]>0)
                {
                    if ( (nr2[j]*peso>E.Grade_Is_Negative_Covered(j) ) or (nr2[j]*peso==E.Grade_Is_Negative_Covered(j) and peso>E.Peso_Is_Negative_Covered(j) ) )
                    {
                        positivos=positivos+(nr2[j]*peso);
                        aciertos++;
                        ejemplos_afectado_por_esta_regla[j]=1;
                        if (nrule>0 and E.Rule_Aislada(j)==-1 and E.Rule_Used(j)>=0)
                          errores_esta_regla[E.Rule_Used(j)]++;
                    }
                    else {
                        positivos_a_incrementar += (nr2[j]*peso);
                    }
                }
                else if (nr2[j]<0) {

                       if (  (-nr2[j]*peso>E.Grade_Is_Positive_Covered(j)) or (-nr2[j]*peso==E.Grade_Is_Positive_Covered(j) and peso>E.Peso_Is_Positive_Covered(j)) ){
                          //negativos=negativos+(-nr2[j]*peso);
                          //if (E.Grade_Is_Positive_Covered(j)>0)
                            negativos=negativos+((-nr2[j]*peso)/ponderacion_por_clase[E.Class_of_example(V.Consecuente(),j)]);

                          ejemplos_afectado_por_esta_regla[j]=-1;

                          if (nrule>0 and E.Rule_Aislada(j)==-1 and E.Rule_Used(j)>=0)
                            errores_esta_regla[E.Rule_Used(j)]++;
                       }

                }
            }
            else
            {
                if (nr2[j]>0)
                {
                   if (nr2[j]*peso>E.Grade_Is_Negative_Covered(j) or
                      (nr2[j]*peso==E.Grade_Is_Negative_Covered(j) and peso>E.Peso_Is_Negative_Covered(j) ) ){
                        ejemplos_afectado_por_esta_regla[j]=1;
                        if (nrule>0 and E.Rule_Aislada(j)==1 and E.Rule_Used(j)>=0)
                           basica_esta_regla[E.Rule_Used(j)]++;
                        positivos_sobrecubrimiento+=1;
                   }

                }
                else  if (nr2[j]<0 and ( (-nr2[j])*peso>E.Grade_Is_Positive_Covered(j) or
                         ((-nr2[j])*peso==E.Grade_Is_Positive_Covered(j) and peso > E.Peso_Is_Positive_Covered(j) ) ) )
                {
                    ejemplos_afectado_por_esta_regla[j]=-1;
                    //negativos=negativos+(-nr2[j]*peso);
                    negativos=negativos+((-nr2[j]*peso)/ponderacion_por_clase[E.Class_of_example(V.Consecuente(),j)]);
                    fallos++;
                    if (nrule>0 and E.Rule_Aislada(j)==1 and E.Rule_Used(j)>=0)
                      basica_esta_regla[E.Rule_Used(j)]++;
                }
            }
        }
    }
// La regla es invalida
    else
    {
        for (j=0; j<ne; j++)
            nr2[j]=0;
        fallos=2*ne;

    }

//delete [] aux;
    //delete [] aux_p;
    //delete [] aux_n;

 if (regla_valida and valida) {

  // Evaluar si es posible eliminar alguna regla ya existente
  int reglas_eliminables=0;
  int reglas_otra_clase=0;
  double aciertos_adicionales=0;
  int ganancia_en_variables=(variables_relevantes+funciones_relevantes);
  int ganancia_funciones_en_variables=funciones_relevantes;
  double diferencia_peso=peso;
  int diferencia_estables=variables_estables+funciones_estables;
  int regla_afectada;
  double max_porcen=0;
  double balance;
  int *lista_reglas_eliminables= new int[nrule];

  for (int i=0; i<nrule; i++){
    if (/*(basica_para_regla[i]==0) or (basica_para_regla[i]-basica_esta_regla[i]==0) or */
       (basica_esta_regla[i]-basica_para_regla[i])+ (errores_para_regla[i]-errores_esta_regla[i])>=0){

      if (true /*reglas_eliminables==0 or R.Class_of_rule(i)!=clase*/) {
        lista_reglas_eliminables[reglas_eliminables]=i;
        reglas_eliminables++;
        if (R.Class_of_rule(i)!=clase)
          reglas_otra_clase++;

        //if (R.Class_of_rule(i)==clase)
          aciertos_adicionales+=(basica_esta_regla[i]-basica_para_regla[i])+ (errores_para_regla[i]-errores_esta_regla[i]);
        /*else
          aciertos_adicionales+=(errores_para_regla[i]-basica_esta_regla[i])-(basica_para_regla[i]-errores_esta_regla[i]);*/

        ganancia_en_variables=R.Variables_of_rule(i)-ganancia_en_variables;
        ganancia_funciones_en_variables=R.New_Atributes_of_rule(i)-ganancia_funciones_en_variables;
        diferencia_peso-=R.Get_Weight(i);
        diferencia_estables-=R.Stable_of_rule(i);
        balance=Diferencia_Adaptacion_Negativa(nr2,R.SeeAdaptaciones_of_rule(i),ne);
        regla_afectada=i;
      }
      else {
        lista_reglas_eliminables[reglas_eliminables]=i;
        reglas_eliminables++;
        int incremento_temporal_aciertos= aciertos_adicionales;
        aciertos_adicionales=(basica_esta_regla[i]-basica_para_regla[i])+ (errores_para_regla[i]-errores_esta_regla[i]);
        double *adapt_nueva_eliminada;
        int j=0;
        while (j<ne and aciertos_adicionales>=0){
          if (ejemplos_afectado_por_esta_regla[j]==0 and E.Rule_Aislada(j)>0 and E.Rule_Aislada(j)<=reglas_eliminables-reglas_otra_clase){
            int k=0, tope=E.Rule_Aislada(j), cont=0;
            while (k<reglas_eliminables and cont+(reglas_eliminables-reglas_otra_clase-k)>=tope){
              if (R.Class_of_rule(lista_reglas_eliminables[k])==clase){
                adapt_nueva_eliminada = R.SeeAdaptaciones_of_rule(lista_reglas_eliminables[k]);
                if (adapt_nueva_eliminada[j]>E.Grade_Is_Negative_Covered(j) or (adapt_nueva_eliminada[j]>0 and adapt_nueva_eliminada[j]==E.Grade_Is_Negative_Covered(j) and R.Get_Weight(lista_reglas_eliminables[k])>=E.Peso_Is_Negative_Covered(j))    ) {
                  cont++;
                }
              }
              k++;
            }
            if (cont==tope)
                aciertos_adicionales--;
          }
          j++;
        }

        if (aciertos_adicionales<0){
          reglas_eliminables--;
          aciertos_adicionales=incremento_temporal_aciertos;
        }
        else {
          aciertos_adicionales+=incremento_temporal_aciertos;
        }

      }

    }
    else  /*if (basica_esta_regla[i]>0 and basica_para_regla[i]-basica_esta_regla[i]!=0)*/
             if (R.Class_of_rule(i)==clase and (1.0*basica_esta_regla[i]+errores_para_regla[i])/(basica_para_regla[i]+errores_esta_regla[i])>max_porcen)
               max_porcen= (1.0*basica_esta_regla[i]+errores_para_regla[i])/(basica_para_regla[i]+errores_esta_regla[i]);
  }

  /*if (reglas_eliminables>=2){
   cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n" << endl;
   for (int i=0; i<nrule; i++){
    if ( (basica_esta_regla[i]-basica_para_regla[i])+ (errores_para_regla[i]-errores_esta_regla[i])>=0 ) ){

      cout << "\t\t\t\tRegla: " << i << endl;
      cout << "\t\t\t\t\tbasica esta regla: "<< basica_esta_regla[i] << endl;
      cout << "\t\t\t\t\tbasica para regla: "<< basica_para_regla[i] << endl;
      cout << "\t\t\t\t\terrores para regla: "<< errores_para_regla[i] << endl;
      cout << "\t\t\t\t\terrores esta regla: "<< errores_esta_regla[i] << endl;
      cout << "\t\t\t\t\tbasica(esta-para): "<< basica_esta_regla[i]-basica_para_regla[i] << endl;
      cout << "\t\t\t\t\terrores(para-esta): "<< errores_para_regla[i]-errores_esta_regla[i] << endl;
      cout << "\t\t\t\t\t\t Aciertos Adicionales: "<< (basica_esta_regla[i]-basica_para_regla[i])+ (errores_para_regla[i]-errores_esta_regla[i]) << endl << endl;
      //sleep(2);
    }
   }
   cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n" << endl;
  }*/

  //double granularidad = 2.0/ponderacion_por_clase[clase];

  double granularidad = 0;
  if (ponderacion_por_clase[clase]>0)
    positivos = positivos/ponderacion_por_clase[clase];

  if (aciertos+aciertos_adicionales-fallos>granularidad){
    if (ponderacion_por_clase[clase]>0) {
      if (positivos-negativos<0)
        valor[0]=-0.9;
      else
        //valor[0]=(positivos/ponderacion_por_clase[clase])-(negativos/(ne-ponderacion_por_clase[clase]));
        //valor[0]=(positivos-negativos)/(ponderacion_por_clase[clase]);
        valor[0]=positivos-negativos;

      if (positivos==0 && menos_infinito)
        valor[0]=-2*ne;
    }
    else
      valor[0]=-2*ne;
  }
  else if (aciertos+aciertos_adicionales>granularidad){
           if (!reducir_reglas){
             valor[0]=aciertos+aciertos_adicionales-fallos-granularidad;
             if (valor[0]==0)
                valor[0]=-ne;
           }
           else {
             valor[0]=aciertos+aciertos_adicionales-fallos-granularidad;
             if (valor[0]==0)
                valor[0]=-ne;
             else if ( reglas_eliminables>=(-valor[0])+2)
                    valor[0]=aciertos+aciertos_adicionales-fallos;
                 else
                   valor[0]=-2*ne;
           }
       }
       else
           //valor[0]=aciertos-fallos;
           valor[0]=-(ne+fallos);


   valor[6]= aciertos_adicionales;
   valor[5]= aciertos-fallos;
   valor[4]= variables_estables/*+(0.00001)*unos_variables*/;
   valor[3]= (variables_relevantes+funciones_relevantes);
   //valor[2]= -negativos_global;
   valor[2]=0;
   valor[1]= reglas_eliminables+(afecta*0.99);
   /*if (reglas_eliminables>=1 and valor[0]>0){
      valor[2]=valor[0];
      valor[0]=valor[0]+((1.0*aciertos_adicionales)/ponderacion_por_clase[clase]);
      valor[4]=reglas_otra_clase;
   }*/

   if ( (reglas_eliminables>=1 and reglas_eliminables<2) and (valor[0]==-ne or valor[0]==0) ){
     if (!discriminante){
       valor[1]=1;
       valor[2]=balance;
       valor[3]=ganancia_en_variables;
       valor[4]=0;
       //valor[4]=ganancia_funciones_en_variables;
       valor[5]=diferencia_estables;
       valor[6]=max_porcen;
     }
     else {
       valor[1]=1;
       //valor[3]=0;
       if (aciertos_adicionales>0)
         valor[2]=0.1+aciertos_adicionales;
       else
         if (balance>0 and ganancia_en_variables>0){
           valor[2]=1;
         }
         else {
           valor[2]=0;
         }
       valor[3]=ganancia_en_variables;
       valor[4]=ganancia_funciones_en_variables;
       valor[5]=diferencia_estables;
       valor[6]=max_porcen;
       //valor[2]=balance;
     }
   }
   delete [] lista_reglas_eliminables;
  }
  else {
    valor[0]=valor[1]=valor[2]=valor[3]=valor[4]=valor[5]=valor[6]=-2*ne;
    peso=0;
  }


  if (nrule>0){
    delete [] basica_esta_regla;
    delete [] errores_esta_regla;
  }

  delete [] ejemplos_afectado_por_esta_regla;
  delete [] nr2;

  return peso;
}
