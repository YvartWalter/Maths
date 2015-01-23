///
/// \version 1.0
/// \date 03/06/2014
/// \author Walter Yvart
///

///
/// Fournir des foncteurs sortran une fonction ou sa dérivée 
///
#ifndef _BIG_B_MATH_MODULE_
static const double Class_Maths_B_var1 = 0.01;
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
    double ix{(double)i + Class_Maths_B_var1};
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
    inline virtual T operator[](unsigned int i) const { if(i<par.size()) return par[i]; else return T{0}; };
    inline virtual const T* Parameters(void) const { return par.data(); };
    inline virtual operator const T* () const { return par.data(); };
    inline virtual operator unsigned int () const { return par.size();};
    inline virtual operator unsigned int () { return par.size();};
    inline virtual T* Parameters(void) { return par.data(); };
    inline virtual operator T* () { return par.data(); };
    
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
/// Fournir une base linéaire automatique ne passant pas par 0 et finissant proche de 1
///
template<typename T>
class auto_base_lin
{
  protected:
    std::vector<T> abs;
    T max;
  public:
    auto_base_lin(void):abs(), max(0) { };
    auto_base_lin(unsigned int len):abs(len), max(T{len})
    {
      while(len!=0){  abs[len-1] = T{Class_Maths_B_var1}+T{len-1}/max; len--; }
    };
    auto_base_lin(const auto_base_lin& a):abs(a.abs), max(a.max) { };
    template<typename U> auto_base_lin(const auto_base_lin<U>& a):abs(), max(0)
    {
      for(unsigned int i=0;i<a.abs.size();++i)
      {
        abs.push_back( T{a.abs[i]});
      }
      max=T(abs.size());
    };
    virtual ~auto_base_lin() { };
    
    void Rebase(unsigned int len)
    {
      abs.clear();
      max=static_cast<T>(len);
      for(unsigned int i=0;i<len;++i)
      {
        abs.push_back( T{Class_Maths_B_var1}+static_cast<T>(i)/max );
      }
    };
    
    void LasyRebase(unsigned int len){  if(len!=abs.size()) Rebase( len );  };
    
    void StrongRebase( T* x, unsigned int len )
    {
      abs.clear();
      max=static_cast<T>(len);
      for(unsigned int i=0;i<len;++i){   abs.push_back( x[i] );   }
    }
    
    inline T at( unsigned int i) const { return abs[i%abs.size()]; };
    inline T at( unsigned int i){ return abs[i%abs.size()]; };
    
    inline T* Base(void) const { return abs.data(); };
    inline T* Base(void){ return abs.data(); };
    
    inline operator unsigned int() { return abs.size(); };
    inline operator unsigned int() const { return abs.size(); };
    
    inline operator bool () {return true; };
    inline operator bool () const {return true; };
};

///
/// Polynome de degré quelconque, le coefficient du degré le plus bas en premier BigEndian
///
template<typename T>
class Polynome : public virtual ParametricBase<T>, private virtual auto_base_lin<T>
{
  public:
    ///
    /// Constructeurs
    ///
    Polynome(void):ParametricBase<T>(), auto_base_lin<T>() { this->par.push_back( T{1} ); };
    Polynome(T* a, unsigned int len):ParametricBase<T>(), auto_base_lin<T>()
    {
      unsigned int u{0};
      while(u!=len)
      {
        if( a!=0) this->par.push_back( a[u] );
        else this->par.push_back( T{1}/T{(u+1)} );
        ++u;
      }
    };
    Polynome(unsigned int len):ParametricBase<T>(), auto_base_lin<T>()
    {
      unsigned int u{0};
      while(u!=len)
      {
        this->par.push_back( T{1}/T{((double)u+1)} );
        ++u;
      }
    };
    Polynome( const Polynome& p ):ParametricBase<T>(p), auto_base_lin<T>(p) { };
    template<typename U> Polynome( const Polynome<U>& p ):ParametricBase<T>(p), auto_base_lin<T>(p) { };
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
    
    void ForAll(T* a, unsigned int len)
    {
      this->auto_base_lin<T>::LasyRebase( len );
      for(unsigned int i=0;i<len;++i)
      {
        a[i] = this( this->auto_base_lin<T>::at(i) );
      }
    };
    
    void ForAll(T* x, T* y, unsigned int len)
    {
      this->auto_base_lin<T>::StrongRebase( x, len );
      for(unsigned int i=0;i<len;++i)
      {
        y[i] = this( x[i] );
      }
    };
    
    ///
    /// Effectuer une inversion en utilisant la fonction LinFit
    ///
    template<class EigenMatrix, class EigenVector>
    T Retrieve(T* x, T* y, T* res, unsigned int len)
    {
      double s;
      T* X{0};
      if(x==0){  this->auto_base_lin<T>::LasyRebase( len ); }
      else{      this->auto_base_lin<T>::StrongRebase( x, len ); }
      X=this->auto_base_lin<T>::Base();
      if( res == 0 )
      {
        s = LinFit<Polynome<T>, EigenMatrix, EigenVector>( *this, len, X, y, this->par.size(),
                                                           this->par.data(), this->par.data() );
      }
      else
      {
        s = LinFit<Polynome<T>, EigenMatrix, EigenVector>( *this, len, X, y, this->par.size(),
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
    
    inline virtual operator unsigned int () const { return this->par.size();};
    inline virtual operator unsigned int () { return this->par.size();};
    
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




