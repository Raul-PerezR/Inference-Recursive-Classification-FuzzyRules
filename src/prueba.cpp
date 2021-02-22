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
#include "patrones.h"

using namespace std;


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


void PintaDiferenciaPatrones(const string &cad1, const string &cad2);
int Comparar2patrones(const string &cad1, const string &cad2);
int Comparar2patronesPinta(const string &cad1, const string &cad2);

string Interseccion(const string &cad1, const string &cad2);
string Union(const string &cad1, const string &cad2);
int Grados_Libertad(string cad);
int DistanciaHamming(const string &s1, const string &s2);
int DiferenciaDistanciaHamming(const string &s1, const string &s2, int pos_inicial);




bool vectoresIguales_E(double *v1, double *v2, int tama);
bool vectoresIguales(double *v1, double *v2, int tama);



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
		return 0;
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
	return peso;
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
		cout << "loading data file ..." << endl;
		while (Existe(fich_extension))
		{
			E.AddExampleFich(fich_extension.c_str(),item_par,x); // Ultimo componente false=orden original, true=orden aleatorio
			cout << "...... loading " << fich_extension << endl;
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

	cout << "Loading process finished ..." << endl;
	fichname = extract_noun(fichname);



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
		if (false) {
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
	}


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
