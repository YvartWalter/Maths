///
/// \version 1.0
/// \date 03/06/2014
/// \author Walter Yvart
///

///
/// Fournir des foncteurs sortran une fonction ou sa dérivée 
///
#ifndef _BIG_B_MATH_MODULE_

///
/// \brief fit une fonction paramétrique F_k(x)=y passée par foncteur, F est linéaire en k
/// \param n_par nombre de paramètres pris pat F
/// \param par valeur de retour des paramètres
/// \param fpar valeur initiale des paramètres
///
template<class Foncteur, class EigenMatrix, class EigenVector>
double LinFit( Foncteur& F, unsigned int n_data, double* x, double* y,
               unsigned int n_par, double* fpar)
{
  return LinFit<Foncteur, EigenMatrix, EigenVector>(F, n_data, x, y, n_par, fpar, 0, 0, 0, true);
}

///
/// \brief Vecrsion avec vecteur F_sol(X)
/// \param par vecteur de retour pour les paramètres restitués, il n'y a pas de risque d'aliasing entre
/// fpar et par, les deux pointeurs peuvent être égaux
/// le paramètre islin est true par défaut
/// 
template<class Foncteur, class EigenMatrix, class EigenVector>
double LinFit( Foncteur& F, unsigned int n_data, double* x, double* y,
               unsigned int n_par, double* fpar, double* par)
{
  return LinFit<Foncteur, EigenMatrix, EigenVector>(F, n_data, x, y, n_par, fpar, par, 0, 0, true);
}

///
/// \brief Version avec plus de paramètres
/// \param islin F est linéaire dans tous ses paramètres, on peut rescaler les abscisses
///
template<class Foncteur, class EigenMatrix, class EigenVector>
double LinFit( Foncteur& F, unsigned int n_data, double* x, double* y,
               unsigned int n_par, double* fpar, double* par,
              double* sigma, double* y_p, bool islin )
{
  EigenMatrix G(n_data, n_par), Cd(n_data, n_data), K,A;
  EigenVector M(n_par), Mp(n_par), D(n_data), Do(n_data), Y, X(n_data), sol(n_par);
  double relative_error, co{1.0}, xmax{1.0};
  //
  Cd.setZero();
  if( y==0 ) return -1.0;
  if(islin)
  {
    if(x!=0){ xmax=x[n_data-1]; }
    else{ xmax = static_cast<double>(n_data); }
  }
  //std::cout<<xmax<<"  "<<islin<<std::endl;
  for(int i=0;i<n_data;++i)
  {
    double ix{(double)i};
    Do(i) = y[i];
    //ABSCISSES
    if(x!=0){     X(i) = x[i] / xmax;     }
    else{         X(i) = ix/xmax;    }
    //SIGMA
    if(sigma==0){   Cd(i,i) = 1.0;    }
    else{           Cd(i,i) = sigma[i]*sigma[i];     }
    //TRANSPOSEE DU JACOBIEN
    //std::cout<<"    ";
    for(int j=0;j<n_par;++j)
    {
      
      if(i==0)
      {
        if(fpar!=0){ M(j) = fpar[j]; }//FG
        else{ M(j) = 0.0; }
      }
      G(i,j) = F.Ppartial( j, X(i) );
      //std::cout<<"    ("<<F[j]<<" "<<X(i)<<")  "<<G(i,j);
    }
    //std::cout<<std::endl;
  }
  //Kt*first guess
  //std::cout<<"fg : \n"<<M<<std::endl;
  D = G*M;
  // Linear least square solving
  // Inversion linéaire, pas d'irétaions
  K = G.transpose()*Cd;  Y = K*(Do-D);  A = K*G;  sol = A.fullPivHouseholderQr().solve(Y);  M = M+sol;
  //Recalcule le F_new(X) correspondant à la solution, utilise x[] au lieu de X()
  F.Redefine(M.data(), M.size());
  //D = G*M;
  for(int i=0;i<D.size();++i)
  {
    D(i) = F(X(i));
  }
  relative_error = D.norm() / Do.norm();
  if(y_p!=0){  for(int i=0;i<n_data;++i){   y_p[i] = D(i);   }  }
  if( par!=0 )
  {
    if(islin)
    {
      co=1.0;
      for(int i=0;i<n_par;++i)
      {
        par[i] = F[i]/co;// = M(i)/co;
        co*=xmax; //remettre les coefficients dans la base des x[] et non plus d X()
      }
    }
    else{  for(int i=0;i<n_par;++i){   par[i] = F[i];   }  }
  }
  return relative_error;
}

// ------------------------------------------------------------------------------------- //
//                                                                                       //
// ------------------------------------------------------------------------------------- //

///
/// Polymorphe, fonctions paramétriques
///
template<class Policy>
class Function : virtual public Policy
{
  public:
    Function(void):Policy() { };
    template<class A> Function(A a1):Policy(a1) { };
    template<class A,class B > Function(A a1, B b1):Policy(a1, b1) { };
    template<class A,class B, class C > Function(A a1, B b1, C c1):Policy(a1, b1, c1) { };
    
    Function(const Function<Policy>& f):Policy(f) { };
    template<class Military> Function(const Function<Military>& f):Policy(f) { };
    
    virtual ~Function() {};
};

template<typename T>
class ParametricBase
{
  protected:
    std::vector< T > par;
  public:
    ///
    /// Constructeurs
    ///
    ParametricBase(void):par() { };
    ParametricBase( unsigned int s ):par(s, T{0}) {  };
    ParametricBase( unsigned int s, const T& a ):par(s, a) {  };
    ParametricBase( int s ):par( abs(s), T{0} ) {  };
    ParametricBase( int s, const T& a ):par( abs(s), a ) {  };
    ParametricBase( const ParametricBase<T>& c ): par(c.par) {  };
    template<typename U> ParametricBase( const ParametricBase<U>& c ): par()
    {  for(auto a=c.par.begin(); a!=c.par.end();++a) par.push_back( T{*a} );  };
    virtual ~ParametricBase() { par.clear(); };
    
    ///
    /// Accéder aux parametres
    ///
    virtual T operator[](unsigned int i) const { if(i<par.size()) return par[i]; else return T{0}; };
    virtual const T* Parameters(void) const { return par.data(); };
    virtual operator const T* () const { return par.data(); };
    virtual operator unsigned int () { return par.size();};
    virtual T* Parameters(void) { return par.data(); };
    virtual operator T* () { return par.data(); };
    
    ///
    /// MLodifier les parametres
    ///
    virtual ParametricBase<T>&  Redefine(T* a, unsigned int len)
    {
      if( len == 0 ) return *this;
      par.clear();
      if( a!=0 )
      {
        for(unsigned int u=0;u<len;++u) par.push_back( a[u] );
      }
      else
      {
        //Comportement par défaut pour les pmolynomes
        for(unsigned int u=0;u<len;++u) par.push_back( T{1}/T{((double)u+1)} );
      }
      return *this;
    };
    
    virtual ParametricBase<T>&  Redefine(std::vector<T>& P) { par = P;};
};

///
/// Polynome de degré quelconque, le coefficient du degré le plus bas en premier BigEndian
///
template<typename T>
class Polynome : public virtual ParametricBase<T>
{
  public:
    ///
    /// Constructeurs
    ///
    Polynome(void):ParametricBase<T>() { this->par.push_back( T{1} ); };
    Polynome(T* a, unsigned int len):ParametricBase<T>()
    {
      unsigned int u{0};
      while(u!=len)
      {
        if( a!=0) this->par.push_back( a[u] );
        else this->par.push_back( T{1}/T{(u+1)} );
        ++u;
      }
    };
    Polynome(unsigned int len):ParametricBase<T>()
    {
      unsigned int u{0};
      while(u!=len)
      {
        this->par.push_back( T{1}/T{((double)u+1)} );
        ++u;
      }
    };
    Polynome( const Polynome& p ):ParametricBase<T>(p) { };
    template<typename U> Polynome( const Polynome<U>& p ):ParametricBase<T>(p) { };
    virtual ~Polynome() { };
    
    ///
    /// Renvoie le résultat de l'opération pour x
    ///
    T operator()(const T& x) const
    {
      T sum{0}, X{1};
      for(unsigned int i=0;i<this->par.size();++i)
      {
        sum = sum + this->par[i] * X;
        X = X * x;
      }
      return sum;
    };
    
    ///
    /// Renvoie la valeur de la dérivée en x
    ///
    T prime(const T& x) const
    {
      T sum{0}, X{1};
      for(unsigned int i=1;i<this->par.size();++i)
      {
        sum = sum + T{i} * this->par[i] * X;
        X = X * x;
      }
      return sum;
    };
    
    ///
    /// Dérivée d-ieme par rapport à 0
    ///
    T derivative(unsigned int d, const T& x) const
    {
      T sum{0}, X{1}, A;
      if(d>=this->par.size()) sum;
      for(unsigned int i=d;i<this->par.size();++i)
      {
        A = T{i};
        for(unsigned int j=1;j<i;++j)
        {
          A=A*T{(i-j)};
        };
        sum = sum + A * this->par[i] * X;
        X = X * x;
      }
      return sum;
    };
    
    ///
    /// dérivée partielle par rapport au j-eme coefficient du polynome en x
    ///
    T Ppartial( unsigned int j, const T& x) const
    {
      if( j>=this->par.size() )return T{0};
      if( j==0 ) return this->par[0];
      if( j==1 ) return x;
      return static_cast<T>( pow(x, j) );
      
    };
    
    ///
    /// Effectuer une inversion en utilisant la fonction LinFit
    ///
    template<class EigenMatrix, class EigenVector>
    T Retrieve(T* x, T* y, T* res, unsigned int len)
    {
      double s;
      if( res == 0 )
      {
        s = LinFit<Polynome<T>, EigenMatrix, EigenVector>( *this, len, x, y, this->par.size(),
                                                           this->par.data(), this->par.data() );
      }
      else
      {
        s = LinFit<Polynome<T>, EigenMatrix, EigenVector>( *this, len, x, y, this->par.size(),
                                                           this->par.data(), this->par.data(),
                                                           0, res, true );
      }
      return T{s};
    };
    
    ///
    /// Modifier les paramètres de la fonction
    ///
    
    ///
    /// Accéder aux paramètres et aux nombres de paramètres
    ///
    Polynome<T>& Parameters(T* a, unsigned int len)
    {
      if( a==0 ) return *this;
      for(unsigned int u=0;u<len;++u) a[u] = this->par[u];
    };
    inline virtual operator bool () const {return true; }; //linéaire en alpha, r = true
    inline virtual operator bool () {return true; }; //linéaire en alpha, r = true
    
};

///
/// y = alpha / x**r
/// linéaire en alpha
/// non linéaire en r
///
template<typename T>
class InversPowerProp : public virtual ParametricBase<T>
{
  public:
    /// Constructeuers
    /// par défaut F : x -> 1/x
    ///
    InversPowerProp(void):ParametricBase<T>() { this->par.push_back( T{1} ); this->par.push_back( T{1} ); };
    InversPowerProp(const T& a, const T& p):ParametricBase<T>() { this->par.push_back( T{1} ); this->par.push_back( T{1} );};
    InversPowerProp(const InversPowerProp<T>& c):ParametricBase<T>(c) { };
    template<typename U> InversPowerProp(const InversPowerProp<U>& c):ParametricBase<T>(c) { };
    virtual ~InversPowerProp() { };
    
    ///
    /// Renvoie le résultat de l'opération pour x, renvoie 0 au point singulier
    ///
    T operator()(const T& x) const
    {
      if(x!=T{0} ) return this->par[0] / pow( x, this->par[1] );
      return T{0};
    };
    
    ///
    /// Renvoie la valeur de la dérivée en x, renvoie 0 au point singulier
    ///
    T prime(const T& x) const
    {
      if(x!=T{0} ) return -(this->par[1]) * this->par[0] / static_cast<T>( pow( x, (this->par[1])+1.0 ) );
      return -T{0};
    };
    
    ///
    /// Calcule la dérivée n-ieme de la fonction en x, retourne 0 au point singulier
    ///
    T derivative(unsigned int d, const T& x) const
    {
      T R, u, v, n{1};
      R = this->par[1];
      u = static_cast<T>( pow( x, this->par[1] )); 
      v = n;
      while( d!=1){   v = v*R;   u = u*x;   R+=n;   --d;   }
      u = u*x;
      if( u!=T{0} )return this->par[0]*v/u;
      return T{0};
    };
    
    ///
    /// Renvoie la dérivée partielle par rapport au paramètre j au point x
    /// Renvoie 0 au point singulier
    ///
    T Ppartial( unsigned int j, const T& x) const
    {
      j=j%2;
      if( x==T{0} ) return T{0};
      if(j==0)
      {
        return T{1} / static_cast<T>( pow( x, this->par[1]));
      }
      else
      {
        if( x > T{0} )
        {
          return static_cast<T>( pow(x, this->par[1]) * log( x ) );
        }
        else
        {
          return T{0}; //root
        }
      }
    };
    
    ///
    /// Fonction d'inversion, donné un nuage de point (x,y), va trouver
    /// les paramètres alpha et r fitant sur le poiint moyen (une solution ...)
    /// modifie les paramètres, donne un bon point de départ
    ///
    T Retrieve(T* x, T* y, unsigned int len)
    {
      unsigned int u{0}, v{len-1}, nu{0}, nv{0};
      T x1, x2, y1, y2;
      T diff, zero={0};
      if(!len || (x==0) || (y==0) ) return zero;
      while( u!=len/2 )
      {
        if( (x[u]>zero) && (y[u]!=zero) ){  x1+=x[u]; y1+=y[u]; ++nu;  }
        if( (x[v]>zero) && (y[v]!=zero) ){  x2+=x[v]; y2+=y[v]; ++nv;  }
        ++u; --v;
      }
      if( !nu || !nv ) return zero;
      x1=x1/T{nu};
      y1=y1/T{nu};
      x2=x2/T{nv};
      y2=y2/T{nv};
      // ON résoud simplement 
      this->par[1] = static_cast<T>( log(y2/y1) / log(x1/x2) );
      this->par[0] = y1 * static_cast<T>( pow(x1, this->par[1]) );
      u=0;
      diff=zero;
      while(u!=len)
      {
        x1 = y[u] - (this->par[0] / pow( x[u], this->par[1] ));
        diff+= static_cast<T>( sqrt(x1*x1) );
        ++u;
      }
      diff=diff/T{len};
      return diff;
    };
    
    ///
    /// Modifier les paramètres de la fonction
    ///
    InversPowerProp<T>& Redefine(const T& a, const T& p){ this->par[0]=a; this->par[1]=p;  return *this; };
    
    ///
    /// Accéder aux paramètres et aux nombres de paramètres
    ///
    InversPowerProp<T>& Parameters(T& a, T& p){  a=this->par[0];  p=this->par[1];  return *this; };
    inline virtual operator bool () const {return false; }; //linéaire en alpha, r = true
    inline virtual operator bool () {return false; }; //linéaire en alpha, r = true
    
};

typedef Function<InversPowerProp<double>> d_invpow;
typedef Function<Polynome<double>> d_polynom;

#define _BIG_B_MATH_MODULE_
#endif

/*
 * Test unitaire
#               d_invpow                                      reference                                                    d_polynom
#1 2    3         4         5             6                7  8    9         10  11        12                   13 14   15        16        17        18
#                                                          i  X(i)  Y(i)  Sol  F_sol(X)(i)  Y0(i)
#  plot "./test.dat" u 2:6 w l lt 2 lw 5,"" u 2:3 w p lc 1, "" u 2:4 w lp lc 1, "" u 2:9 w p lc -1, "" u 2:11 w l lt -1 lw 2, "" u 2:15 w p lc 3, "" u 2:16 w l lt 1 lc 3 lw 2
#
0  0.1  139.287   556.878   6.39678e+10   147.419          0  0.1  137.411   0   146.063   147.419              0  0.1  138.419   138.002   138.002   147.419
1  0.2  122.72   291.406   9.30246e+07   134.429           1  0.2  143.708   0   135.652   134.429              1  0.2  134.787   129.593   129.593   134.429
2  0.3  121.445   199.512   2.03632e+06   122.749          2  0.3  131.659   0   125.811   122.749              2  0.3  116.161   121.504   121.504   122.749
3  0.4  116.211   152.488   135280   112.928               3  0.4  120.487   0   116.526   112.928              3  0.4  101.802   113.737   113.737   112.928
4  0.5  99.5047   123.791   16512.3   104.412              4  0.5  114.707   0   107.784   104.412              4  0.5  108.816   106.291   106.291   104.412
5  0.6  105.844   104.402   2961.31   96.7632              5  0.6  92.3243   0   99.5706   96.7632              5  0.6  99.8007   99.1662   99.1662   96.7632
6  0.7  90.1807   90.3976   692.588   89.8072              6  0.7  85.6096   0   91.8725   89.8072              6  0.7  92.8095   92.3626   92.3626   89.8072
7  0.8  75.5685   79.7946   196.73   83.4222               7  0.8  76.3214   0   84.6761   83.4222              7  0.8  89.6948   85.8801   85.8801   83.4222
8  0.9  75.7405   71.4792   64.8234   77.5066              8  0.9  61.4193   0   77.9679   77.5066              8  0.9  83.1395   79.7188   79.7188   77.5066
9  1  77.1007   64.7779   24.0129   71.9985                9  1  71.4153   0   71.7343   71.9985                9  1  80.2823   73.8787   73.8787   71.9985
10  1.1  63.5262   59.2587   9.77904   66.8701             10  1.1  72.469   0   65.9618   66.8701              10  1.1  76.6542   68.3597   68.3597   66.8701
11  1.2  75.4157   54.6318   4.30645   62.1014             11  1.2  77.3232   0   60.6366   62.1014             11  1.2  62.024   63.1619   63.1619   62.1014
12  1.3  58.7042   50.6951   2.0252   57.6506              12  1.3  55.6751   0   55.7452   57.6506             12  1.3  64.9502   58.2853   58.2853   57.6506
13  1.4  40.783   47.3037   1.00719   53.4605              13  1.4  60.929   0   51.2741   53.4605              13  1.4  38.4685   53.7298   53.7298   53.4605
14  1.5  50.3715   44.3506   0.525645   49.4906            14  1.5  45.0774   0   47.2096   49.4906             14  1.5  47.0827   49.4955   49.4955   49.4906
15  1.6  47.6419   41.7552   0.286093   45.7377            15  1.6  43.1828   0   43.5383   45.7377             15  1.6  44.755   45.5824   45.5824   45.7377
16  1.7  54.579   39.4558   0.161564   42.2168             16  1.7  38.7046   0   40.2463   42.2168             16  1.7  38.1391   41.9904   41.9904   42.2168
17  1.8  39.6351   37.4039   0.0942689   38.9176           17  1.8  40.1974   0   37.3203   38.9176             17  1.8  46.5171   38.7196   38.7196   38.9176
18  1.9  33.9574   35.5613   0.05663   35.8143             18  1.9  24.307   0   34.7466   35.8143              18  1.9  31.7073   35.77   35.77   35.8143
19  2  33.5271   33.8973   0.0349205   32.8988             19  2  44.6193   0   32.5116   32.8988               19  2  19.9894   33.1416   33.1416   32.8988
20  2.1  26.7121   32.3867   0.0220475   30.1586           20  2.1  19.3861   0   30.6017   30.1586             20  2.1  39.7647   30.8343   30.8343   30.1586
21  2.2  28.367   31.0092   0.0142211   27.5675            21  2.2  20.9497   0   29.0034   27.5675             21  2.2  33.581   28.8482   28.8482   27.5675
22  2.3  21.4248   29.7476   0.00935345   25.1187          22  2.3  31.0969   0   27.703   25.1187              22  2.3  20.9817   27.1833   27.1833   25.1187
23  2.4  25.3911   28.5879   0.00626262   22.9306          23  2.4  28.4544   0   26.687   22.9306              23  2.4  25.3165   25.8395   25.8395   22.9306
24  2.5  14.114   27.5181   0.00426239   21.4945           24  2.5  29.3526   0   25.9418   21.4945             24  2.5  30.2494   24.8169   24.8169   21.4945
*/


