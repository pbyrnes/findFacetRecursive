#include<iostream>
#include<fstream>
#include<ctime>
#include<string>
#include<cstdlib>
#include<cstring>

using namespace std;

//const int cutoff = 50;
//const int depthCutoff = 1;
static int numGetFacetCalls =0;
static int maxDepth =0;
int d;
int numVertices;
static int fileCount = 0;
void getVertices(int, int&, int&);
int getIndex(int a, int b);
int rank(int** points, int numPoints);
int rank2(int** points, int numPoints);
int groupCutoff;
void getMultipliers(int a, int b, int& tmp1, int& tmp2);
static double totSystemTime = 0;

struct Facet
{
  int* a;
  int b;
  Facet* next;
};

struct pointList
{
  int* a;
  pointList* next;
};

struct pointListList
{
  pointList* pL;
  pointListList* next;
  int numPoints;
  Facet fac;
};

pointListList* pLL = NULL;

struct node
{
  int depth;
  node* left;
  node* right;
  int count;
};

bool lexSmaller(int* f, int* g);
int* applyPermutation(int* sigma, int* f);
void addPoint(node* n, int* p);
void addPoint(node* n, int* p, int k);
bool contains(node* n, int* p);
bool contains(node* n, int* p, int k);
void deleteNode(node* n);

struct Equality
{
  int* a;
  int b;
  Equality* next;
};

struct Partition
{
  int* parts;
};

int* rotate(int* g, int** points, int* f, int N);
void lexMinRecursion(int k, int* sigma, int* S, int* sigmaStar, int* f, Partition vPart);
Equality* equalityCopy(Equality* e);
void addEquality(Equality* e, int* a, int b);
Partition refinePartition(Partition part, int** points, int numPoints);
char* runPorta(int** points, int numPoints);
Facet getPortaResults(char* fileName);
int** getPoints(Facet f, int** points, int oldNumPoints, int& newNumPoints);
bool addFacet(Facet& f, Facet* newFacet);
Facet getFacets(int** points, int numPoints, Partition vPart, Equality* 
hyperP,int cutoff);
Facet getFacetsByPorta(int** points, int numPoints, Partition vPart, Equality* hyperP);
int** readPoints(char* fileName, int& N);
void outputFacet(Facet* fPtr);
Partition singleBlockPartition(int n);
Facet getSingleFacet(int** points, int numPoints);
void normalizeFacet(Facet& f, Partition vPart, Equality* hyperP);
int quotientGroupSize(int* parts);
int* getIndices(int** m, int nE);

int main(int argc, char** argv)
{
  time_t start,end;
  time(&start);
  //input points
  int cutoff = atoi(argv[2]);
  groupCutoff = atoi(argv[3]);
  int N;
  int** points = readPoints(argv[1], N);
//  int r = rank(points,N);
  
  Partition initPart = singleBlockPartition(numVertices);
  Equality hypP;
  hypP.a = new int[d];
  for(int i=0; i<d; i++)
    hypP.a[i] = 1;
  hypP.b = numVertices - 1;
  hypP.next = NULL;

  Facet theFacets = getFacets(points, N, initPart, &hypP,cutoff);

  Facet* fPtr = &theFacets;

  cout << endl;
  cout << "Facets found: " << endl;
  while(fPtr != NULL)
  {
    outputFacet(fPtr);
    Facet* oldF = fPtr;
    fPtr = fPtr->next;
    delete[] oldF->a;
    if(oldF != &theFacets)
      delete oldF;
  }
  for(int i=0; i<N; i++)
    delete[] points[i];
  delete[] points;
  delete[] hypP.a;
  delete[] initPart.parts;

  time(&end);
  cout << "Total of " << numGetFacetCalls << " calls of getFacets" << 
endl;
  cout << "Maximum depth of : " << maxDepth << endl;
  cout << "Computation took: " << difftime(end,start) << " seconds" << 
endl;
  cout << "System call time totalled: " << totSystemTime << " seconds" << endl;
  return 0;
}

Partition singleBlockPartition(int n)
{
  Partition p;
  p.parts = new int[n+1];
  for(int i=0; i<=n; i++)
  {
    p.parts[i] = 1;
  }
  return p;
}

void outputFacet(Facet* fPtr)
{
  for(int i=0; i<d; i++)
    cout << fPtr->a[i] << " ";
  cout << fPtr->b << endl;
}

/*Facet newFacet(int* f)
{
  Facet newF;
  
  newF.next = NULL;
  newF.a = new int[d];
  for(int i=0; i<d; i++)
    newF.a[i] = f[i];
  newF.b = f[d];
  return newF;
}*/

int** readPoints(char* fileName, int& N)
{
  ifstream fin(fileName);
  if(!fin)
  {
    cout << "Error(1): trouble opening " << fileName << endl;
    exit(1);
  }
  fin >> N;
  fin >> d;
  fin >> numVertices;

  int** p = new int*[N];
  for(int i=0; i<N; i++)
  {
    int* tmp = new int[d];
    p[i] = tmp;
    for(int j=0; j<d; j++)
      fin >> p[i][j];
  }
  fin.close();

  return p;
}

Facet getFacets(int** points, int numPoints, Partition vPart, Equality* 
hyperP, int cutoff)
{
  numGetFacetCalls++;
  if(numGetFacetCalls % 1000 == 0)
    cout << "starting getFacets call number: " << numGetFacetCalls << endl;
  int depth = 0;
  Equality* dPtr = hyperP;
  while(dPtr->next != NULL)
  {
    dPtr = dPtr->next;
    depth++;
  }
  if(depth > maxDepth)
    maxDepth = depth;
  int* passParts = new int[d+1];
  for(int i=0; i<=d; i++)
    passParts[i] = vPart.parts[i];
  int q = quotientGroupSize(passParts);
  delete[] passParts;
  if(numPoints < cutoff || q > groupCutoff)
    return getFacetsByPorta(points, numPoints, vPart, hyperP);

  Facet allFacets = getSingleFacet(points, numPoints);
  normalizeFacet(allFacets, vPart, hyperP);
  if(depth == 0)
    outputFacet(&allFacets);

  Facet* currentFacet = &allFacets;

  while(currentFacet != NULL)
  {
    int newNumPoints;
    int** newPoints = getPoints(*currentFacet,points,numPoints,newNumPoints);
    if(newNumPoints < 2)
    {
      cout << "AAA" << endl;
    }
    Equality* newHyperP;
    newHyperP = equalityCopy(hyperP);
    addEquality(newHyperP,currentFacet->a, currentFacet->b);
    Partition newVPart = refinePartition(vPart, newPoints,newNumPoints);
    int numZeros = 0;
    for(int i=0; i<d; i++)
    {
      if(currentFacet->a[i] == 0)
        numZeros++;
    }
    if(numZeros != d-1 || depth > 0)
    {
      // we don't have an x_i >= 0 facet
      
    Facet ridges = getFacets(newPoints, newNumPoints, newVPart, 
newHyperP,cutoff);
    delete[] newVPart.parts;
    while(newHyperP != NULL)
    {
      delete[] newHyperP->a;
      Equality* ePtr = newHyperP;
      newHyperP = newHyperP->next;
      delete ePtr;
    }
    for(int i=0; i<newNumPoints; i++)
      delete[] newPoints[i];
    delete[] newPoints;
    Facet* fPtr = &ridges;
    while(fPtr != NULL)
    {
      int* currentF = new int[d+1];
      for(int i=0; i<d; i++)
        currentF[i] = currentFacet->a[i];
      currentF[d] = currentFacet->b;
      int* g = new int[d+1];
      for(int i=0; i<d; i++)
        g[i] = fPtr->a[i];
      g[d] = fPtr->b;
      
      int* h2 = rotate(g, points, currentF,numPoints);
      for(int i=0; i<d; i++)
        fPtr->a[i] = h2[i];
      fPtr->b = h2[d];
      int newNumPoints2;
      normalizeFacet(*fPtr,vPart,hyperP);

      delete[] g;
      delete[] h2;
      delete[] currentF;
/*cout << "about to add new facet: ";
outputFacet(fPtr);
cout << endl;*/
      bool newFlag = addFacet(allFacets, fPtr);
      if(depth < 1 && newFlag)
        outputFacet(fPtr);
      Facet* oldFPtr = fPtr;
      fPtr = fPtr->next;
      delete[] oldFPtr->a;
      if(oldFPtr != &ridges)
        delete oldFPtr;
    }
    }
    currentFacet = currentFacet->next;
//    if(depth < 1)
    {
      int numFacetsLeft = 0;
      Facet* fCtr = currentFacet;
      while(fCtr != NULL)
      {
        numFacetsLeft++;
        fCtr = fCtr->next;
      }
      cout << numFacetsLeft << " facets left to process at depth " << depth << endl;
    }
  }

  return allFacets;
}

Equality* equalityCopy(Equality* e)
{
  Equality* newE = new Equality;
  newE->a = new int[d];
  for(int i=0; i<d; i++)
    newE->a[i] = e->a[i];
  newE->b = e->b;
  newE->next = NULL;
  Equality* currEOld = e;
  Equality* currENew = newE;

  while(currEOld->next != NULL)
  {
    currEOld = currEOld->next;
    currENew->next = new Equality;
    currENew = currENew->next;
    currENew->a = new int[d];
    for(int i=0; i<d; i++)
      currENew->a[i] = currEOld->a[i];
    currENew->b = currEOld->b;
    currENew->next = NULL;
  }
  return newE;
}

void addEquality(Equality* e, int* a, int b)
{
  Equality* ePtr = e;
  while(ePtr->next != NULL)
  {
    ePtr = ePtr->next;
  }
  ePtr->next = new Equality;
  (ePtr->next)->a = new int[d];
  (ePtr->next)->next = NULL;
  (ePtr->next)->b = b;
  for(int i=0; i<d; i++)
    (ePtr->next)->a[i] = a[i];
}

Partition refinePartition(Partition p, int** points, int numPoints)
{
  Partition newP;
  node* n = new node;
  n->count = 0;
  n->left = NULL;
  n->right = NULL;
  n->depth = 0;
  for(int i=0; i<numPoints; i++)
  {
    addPoint(n,points[i]);
  }

  int* sigma = new int[numVertices];
  bool** noncompatible = new bool*[numVertices+1];
  for(int i=1; i<=numVertices; i++)
  {
    noncompatible[i] = new bool[numVertices+1];
    for(int j=1; j<=numVertices; j++)
      noncompatible[i][j] = false;
  }
  for(int s1=1; s1<numVertices; s1++)
  {
  for(int s2=s1+1; s2<=numVertices; s2++)
  {

    for(int i=0; i<numVertices; i++)
      sigma[i] = i+1;
    sigma[s1-1] = s2;
    sigma[s2-1] = s1;

  bool flag2 = true;
  
  for(int i=0; i<numPoints && flag2; i++)
  {
      int* tmp = applyPermutation(sigma, points[i]);
      if(!contains(n,tmp))
      {
        noncompatible[s1][s2] = true;
        noncompatible[s2][s1] = true;      
        flag2 = false;
      }
      delete[] tmp;
  }
  }
  }
  
  deleteNode(n);

  newP.parts = new int[numVertices+1];
  for(int i=0; i<=numVertices; i++)
    newP.parts[i] = -1;
  int nextPart = 1;
  for(int i=1; i<=numVertices; i++)
  {
    if(newP.parts[i] == -1)
    {
      newP.parts[i] = nextPart;
      for(int j=i+1; j<=numVertices; j++)
      {
        if(p.parts[i] == p.parts[j] && !noncompatible[i][j])
        {
          if(newP.parts[j] != -1)
            cout << "Error: assigning a vertex to 2 parts of a partition" << endl;
          newP.parts[j] = nextPart;
        }
      }
      nextPart++;
    }
  }
  delete[] sigma;
  for(int i=1; i<=numVertices; i++)
    delete [] noncompatible[i];
  delete[] noncompatible;
  return newP;
}

int* applyPermutation(int* sigma, int* f)
{
  int* tmp = new int[d+1];
  for(int i=0; i<d; i++)
  {
    int a,b;
    getVertices(i,a,b);
    tmp[i] = f[getIndex(sigma[a-1],sigma[b-1])];
  }
  tmp[d] = f[d];
  return tmp;
}

int getIndex(int a, int b)
{
  if(a == b)
    cout << "ERROR(5): a==b in getIndex" << endl;
  int t = 0;
  if(a>b)
  {
    int tmp = a;
    a = b;
    b = tmp;
  }
  for(int i=1; i<a; i++)
  {
    t += numVertices - i;
  }
  t += b-a;
  return t-1;
}

void getVertices(int ind, int& a, int& b)
{
  int t = ind+1;
  a = 1;
  b = 0;
  while(t > numVertices-a)
  {
    t -= numVertices-a;
    a++;
  }
  b = a + t;
}

void deleteNode(node* n)
{
  if(n->left != NULL)
    deleteNode(n->left);
  if(n->right != NULL)
    deleteNode(n->right);
  delete n;
  return;
}

void addPoint(node* n, int* p)
{
  addPoint(n, p, 0);
}

void addPoint(node* n, int* p, int k)
{
  n->count = n->count + 1;
  if(k >= d)
    return;
  if(p[0] == 0)
  {
    if(n->left == NULL)
    {
      node* n2 = new node;
      n->left = n2;
      n2->depth = k+1;
      n2->left = NULL;
      n2->right = NULL;
      n2->count = 0;
    }
    addPoint(n->left,&(p[1]),k+1);
  }
  if(p[0] == 1)
  {
    if(n->right == NULL)
    {
      node* n2 = new node;
      n->right = n2;
      n2->depth = k+1;
      n2->left = NULL;
      n2->right = NULL;
      n2->count = 0;
    }

    addPoint(n->right,&(p[1]),k+1);

  }
  return;
}

bool contains(node* n, int* p)
{
  return contains(n,p,0);
}

bool contains(node* n, int* p, int k)
{
  if(k==d)
  {
    if(n->count > 0)
      return true;
    else
      return false;
  }
  if(p[0] == 0)
  {
    if(n->left == NULL)
      return false;
    else
      return contains(n->left,&(p[1]),k+1);
  }
  if(p[0] == 1)
  {
    if(n->right == NULL)
      return false;
    else
      return contains(n->right,&(p[1]),k+1);
  }
  cout << "ERROR p[0] is not 0 or 1, p[0] is: " << p[0] << endl;
  cout << "k is: " << k << endl;
  return false;
}


char* runPorta(int** points, int numPoints)
{
/*  srand(time(NULL));
  int rnd = rand()%10000000;*/
  char fileName[100];
  sprintf(fileName, "tmp%d.poi",fileCount);

  
  ofstream fout(fileName);
  if(!fout)
  {
    cout << "Error: could not open " << fileName << endl;
    exit(1);
  }
  fout << "DIM = " << d << endl;
  fout << "CONV_SECTION" << endl;
  for(int i=0; i<numPoints; i++)
  {
    for(int j=0; j<d; j++)
    {
      fout << points[i][j] << " ";
    }
    fout << endl;
  }
  fout << "END" << endl;
  fout.close();

  char cmd[100];
  sprintf(cmd,"./traf -p %s", fileName);
  time_t start, end;
  time(&start);
  system(cmd);
  time(&end);
  totSystemTime += difftime(end,start);

  char* returnString = new char[100];
  strcpy(returnString, fileName);
  return returnString;
}

Facet getPortaResults(char* fileName)
{
  Facet f;
  f.a = new int[d];
  Facet* fPtr = &f;
  char inFileName[100];
  sprintf(inFileName,"%s.ieq",fileName);
  ifstream fin(inFileName);
  if(!fin)
  {
    cout << "Error: could not open " << inFileName << endl;
    exit(1);
  }

    const int NN = 200;
    char s[NN];
    s[0] = 0;
    while(strchr(s,'<') == NULL)
    {
      fin.getline(s,NN);
    }
    bool flag2 = true;
    while(flag2)
    {
      int* g = new int[d+1];
      for(int i=0; i<=d; i++)
        g[i] = 0;
      char* ss = strtok(s,"x");
      char* tmpS = strchr(s,')');
      ss = tmpS + sizeof(char);
      int lastCo = atoi(ss);
      if(lastCo == 0)
      {
        //ss is "-" or "+"
        bool flag = true;
        for(int i=0; i<strlen(ss) && flag; i++)
        {
          if(ss[i] == '-')
          {
            flag = false;
            lastCo = -1;
          }
          if(ss[i] == '+')
          {
            flag = false;
            lastCo = 1;
          }
        }
      }
      while(ss = strtok(NULL,"x"))
      {
        int t=0;
        while(ss[t] != '-' && ss[t] != '+' && ss[t] != 0)
          t++;
        char sss[NN];
        for(int i=0; i<t; i++)
          sss[i] = ss[i];
        sss[t] = 0;
        int ind = atoi(sss);
        g[ind-1] = lastCo;
        lastCo = atoi(ss+t*sizeof(char));
        if(lastCo == 0)
        {
          //ss is "-" or "+"
          bool flag = true;
          for(int i=t; i<strlen(ss) && flag; i++)
          {
            if(ss[i] == '-')
            {
              flag = false;
              lastCo = -1;
            }
            if(ss[i] == '+')
            {
              flag = false;
              lastCo = 1;
            }
          }
        }
      }
      int tt = 0;
      for(int i=0; i<NN; i++)
      {
        if(s[i] == '=')
        {
          tt = i+1;
          break;
        }
      }
      g[d] = atoi(s+tt*sizeof(char));
      for(int i=0; i<d; i++)
        fPtr->a[i] = g[i];
      fPtr->b = g[d];
      fPtr->next = new Facet;
      fPtr = fPtr->next;
      fPtr->a = new int[d];
      fPtr->next = NULL;
      delete[] g;
      char c;
      fin.get(c);
      if(c != '(')
      {
        flag2 = false;
      }
      else
      {
        fin.getline(s,NN);
      }
    }

   return f;
}

Facet getSingleFacet(int** points, int numPoints)
{
  int smallNumPoints = 50;
  if(smallNumPoints > numPoints)
    smallNumPoints = numPoints;
  Facet resultF;
  resultF.a = new int[d];
  bool noFacetFound = true;
  int initRank = rank(points,numPoints);
  while(smallNumPoints <= numPoints && noFacetFound)
  {
    char* fileName = runPorta(points,smallNumPoints);
    Facet f = getPortaResults(fileName); 
    Facet* initF = &f;
    delete[] fileName;
    while(f.next != NULL && noFacetFound)
    {
      //check if f is a facet for the whole polytope
      bool stillValid = true;
      int numFPoints = 0;
      int** fPoints = new int*[numPoints];
      for(int i=0; i<numPoints; i++)
        fPoints[i] = new int[d];
      for(int i=0; i<numPoints && stillValid; i++)
      {
        int fp = 0;
        for(int j=0; j<d; j++)
          fp += f.a[j]*points[i][j];
        if(fp > f.b)
          stillValid = false;
        if(fp == f.b)
        {
          for(int j=0; j<d; j++)
            fPoints[numFPoints][j] = points[i][j];
          numFPoints++;
        }         
      }
      int fRank = rank(fPoints, numFPoints);
      if(fRank < initRank-1)
        stillValid = false;
      if(fRank > initRank)
      {
        cout << "Error: rank went up somehow" << endl;
        exit(1);
      }
      for(int i=0; i<numPoints; i++)
        delete[] fPoints[i];
      delete[] fPoints;
      if(stillValid)
      {
        if(numFPoints < 2)
        {
          cout << "BBB" << endl;
        }
        noFacetFound = false;
        for(int i=0; i<d; i++)
          resultF.a[i] = f.a[i];
        resultF.b = f.b;
        resultF.next = NULL;
      }
      else
        f = *(f.next);
    }
    while(initF != NULL)
    {
      delete[] initF->a;
      initF = initF->next;
    }
    if(noFacetFound)
    {
      smallNumPoints = smallNumPoints + 50;
      if(smallNumPoints > numPoints)
        smallNumPoints = numPoints;
    }
  }
  if(noFacetFound)
  {
    cout << "Error: did not find a single facet" << endl;
    exit(1);
  }

  return resultF;
}

void normalizeFacet(Facet& fac, Partition vPart, Equality* hyperP)
{
  int* vertSet = new int[d+1];
  for(int i=1; i<=numVertices; i++)
    vertSet[i] = vPart.parts[i];
  int* f = new int[d+1];
  f[d] = fac.b;
  for(int i=0; i<d; i++)
    f[i] = (fac.a)[i];
  if(f[d] != 0)
  {
    int tmp = f[d];
    for(int i=0; i<d; i++)
    {
      f[i] = f[i] * (numVertices-1)-tmp;
    }
    f[d] = 0; 
  }    
  bool fflag = false;
  int minCo = abs(f[0]);
  for(int i=1; i<=d; i++)
  {
    if(abs(f[i]) > 0)
    {  
      if(!fflag)
      {
        minCo = abs(f[i]);
        fflag = true;
      }
      if(abs(f[i]) < minCo)
        minCo = abs(f[i]);
    }
  }
  bool flag1 = true;
  while(flag1)
  {    
    flag1 = false;
    for(int j=2; j<=minCo; j++)
    {
      bool flag = true;        
      for(int i=0; i<=d && flag; i++)
      {
        if((f[i] % j) != 0)
          flag = false;
      }
      if(flag)
      {
        flag1 = true;
        for(int i=0; i<=d; i++)
        {
          f[i] = f[i] / j;
        }
      }
    }
  }
 
//now f has RHS = 0, all coefficients integer, and gcd(coefficients) = 1
  int numEq = 0;
  Equality* ePtr = hyperP;
  while(ePtr->next !=NULL)
  {
    ePtr = ePtr->next;
    numEq++;
  }
  int** e = new int*[numEq];
  for(int i=0; i<numEq; i++)
    e[i] = new int[d];
  ePtr = hyperP->next;
  for(int i=0; i<numEq; i++)
  {
    for(int j=0; j<d; j++)
      e[i][j] = ePtr->a[j];
    ePtr = ePtr->next;
  }
  int* jA = getIndices(e,numEq);
  for(int k=0; k<numEq; k++)
  {
//    ePtr = ePtr->next;
    int j = jA[k];
  //    if((ePtr->a)[j] != 0)
      {
        if(f[j] != 0)
        {
          int tmp1;
          int tmp2;
          getMultipliers(e[k][j], f[j],tmp1,tmp2);
          if(e[k][j] == 0)
          {
            cout << "Error: invalid equality index" << endl;
            exit(1);
          }
          if(e[k][j] > 0)
          {
            for(int i=0; i<d; i++)
            {
              f[i] = f[i]*tmp1-tmp2*e[k][i];
            }
          }
          else
          {
            for(int i=0; i<d; i++)
            {
              f[i] = -f[i]*tmp1+tmp2*e[k][i];
            }
          }
        }
      }
  }
  delete[] jA;
  for(int i=0; i<numEq; i++)
    delete[] e[i];
  delete[] e;
    
/*
  bool* alreadyZeroed = new bool[d];
  for(int i=0; i<d; i++)
    alreadyZeroed[i] = false;
  Equality* ePtr = hyperP;
  bool zeroFlag;
  while(ePtr->next != NULL)
  {
    ePtr = ePtr->next;
    zeroFlag = true;
    for(int j=0; j<d && zeroFlag; j++)
    {
      if((ePtr->a)[j] != 0 && alreadyZeroed[j] == false)
      {
        alreadyZeroed[j] = true;
        zeroFlag = false;
        if(f[j] != 0)
        {
          if((ePtr->a)[j] > 0)
          {
            int tmp = f[j];
            for(int i=0; i<d; i++)
            {
              f[i] = f[i]*(ePtr->a)[j]-tmp*(ePtr->a)[i];
            }
          }
          else
          {
            int tmp = f[j];
            for(int i=0; i<d; i++)
            {
              f[i] = -f[i]*(ePtr->a)[j]+tmp*(ePtr->a)[i];
            }
          }
        }
      }
    }
  }
    if(zeroFlag)
    {
      ePtr = hyperP;
      while(ePtr->next != NULL)
      {
        ePtr = ePtr->next;
        for(int j=d-1; j>=0 && zeroFlag; j--)
        {
          if((ePtr->a)[j] != 0 && alreadyZeroed[j] == false)
          {
            alreadyZeroed[j] = true;
            zeroFlag = false;
            if(f[j] != 0)
            {
              if((ePtr->a)[j] > 0)
              {
                int tmp = f[j];
                for(int i=0; i<d; i++)
                {
                  f[i] = f[i]*(ePtr->a)[j]-tmp*(ePtr->a)[i];
                }
              }
              else
              {
                int tmp = f[j];
                for(int i=0; i<d; i++)
                {
                  f[i] = -f[i]*(ePtr->a)[j]+tmp*(ePtr->a)[i];
                }
              }
            }
          }
        }
      }
      if(zeroFlag)
      {
        cout << "Error(zF)" << endl;
        for(int i=0; i<=d; i++)
          cout << f[i] << ' ';
        cout << endl;
        exit(1);
      }
    }
  delete[] alreadyZeroed;*/
  fflag = false;
  minCo = abs(f[0]);
  for(int i=1; i<=d; i++)
  {
    if(abs(f[i]) > 0)
    {
      if(!fflag)
      {  
        minCo = abs(f[i]);
        fflag = true;
      }
      if(abs(f[i]) < minCo)
        minCo = abs(f[i]);
    }
  }
  flag1 = true;
  while(flag1)
  {
    flag1 = false;
    for(int j=2; j<=minCo; j++)
    {
      bool flag = true;
      for(int i=0; i<=d && flag; i++)
      {
        if((f[i] % j) != 0)
          flag = false;
      }
      if(flag)
      {
        flag1 = true;
        for(int i=0; i<=d; i++)
        {
          f[i] = f[i] / j;
        }  
      }
    }
  }
  int* sigmaStar;
  int* sigma;
  sigmaStar = new int[numVertices];
  for(int i=0; i<numVertices; i++)
  {
    sigmaStar[i] = i+1;
  }
             
  int min = f[0];
  for(int i=1; i<d; i++)
  {
    //check if we can swap x_0 and x_i
    int a0,b0,ai,bi;
    a0=1;
    b0=2;
    getVertices(i,ai,bi);
    if((vertSet[ai] == vertSet[a0] && vertSet[bi] == vertSet[b0]) ||
       (vertSet[ai] == vertSet[b0] && vertSet[bi] == vertSet[a0]))
    {
      if(f[i] < min) 
      {
        min = f[i];
      }  
    }
  }
       
  for(int i=0; i<d; i++)   
  {
    int a0,b0,ai,bi;
    a0=1;
    b0=2;
    getVertices(i,ai,bi);
    if(f[i] == min && (vertSet[ai] == vertSet[a0] && vertSet[bi] == vertSet[b0]))
    {
      sigma = new int[numVertices];
      for(int j=0; j<numVertices; j++)
        sigma[j] = 0;
      int a,b;
      getVertices(i,a,b);  
      sigma[0] = a;
      sigma[1] = b;
      int* S; 
      S = new int[numVertices];
      for(int j=0; j<numVertices; j++)
        S[j] = 0;
      S[a-1] = 1;
      S[b-1] = 1;
      lexMinRecursion(3,sigma,S,sigmaStar,f,vPart);
      delete[] S;
      delete[] sigma;
    }
   
    if(f[i] == min && (vertSet[ai] == vertSet[b0] && vertSet[bi] == vertSet[a0]))
    {
      sigma = new int[numVertices];
      for(int j=0; j<numVertices; j++)
        sigma[j] = 0;
      int a,b;
      getVertices(i,a,b);
      sigma[0] = b;
      sigma[1] = a;
      int* S;
      S = new int[numVertices];
      for(int j=0; j<numVertices; j++)
        S[j] = 0;
      S[a-1] = 1;
      S[b-1] = 1;
      lexMinRecursion(3,sigma,S,sigmaStar,f,vPart);
      delete[] S;   
      delete[] sigma;
    }
  }
  int* newF =  applyPermutation(sigmaStar,f);
  delete[] sigmaStar;
  delete[] f;
  for(int i=0; i<d; i++)
  {
    fac.a[i] = newF[i];
  }
  fac.b = newF[d];
  delete[] newF;
  delete[] vertSet;
  return;
}

bool lexSmaller(int* f, int* g)
{
  for(int i=0; i<d; i++)
  {
    if(f[i] < g[i])
      return true;
    if(f[i] > g[i])
      return false;
  }
  return false;
}


void lexMinRecursion(int k, int* sigma, int* S, int* sigmaStar, int* f, Partition vPart)
{
  int* vertSet = new int[numVertices+1];
  for(int i=1; i<=numVertices; i++)
    vertSet[i] = vPart.parts[i];
      
  if(k == numVertices)
  {
    bool viableFlag = false;
    for(int i=0; i<numVertices && !viableFlag; i++)
    {
      if(S[i] == 0)
      {
        sigma[k-1] = i+1;
        if(vertSet[k] == vertSet[i+1])
          viableFlag = true;
      }
    }
    int* x = applyPermutation(sigma,f);
    int* y = applyPermutation(sigmaStar,f);
    if(lexSmaller(x, y) && viableFlag)
    {
      bool checkFlag = true;
      for(int i=0; i<numVertices; i++)
      {
        if(vertSet[sigma[i]] != vertSet[i+1])
          cout << "Error: invalid sigmaStar" << endl;
      }
      for(int i=0; i<numVertices; i++)
        sigmaStar[i] = sigma[i];
    }
    delete[] x;
    delete[] y;
  }
  else
  {
    if(k > numVertices)
    {
      cout << "ERROR(4): k > numVertices in call to lexMinRecursion" << endl;
      return;
    }
    bool flag = false;
    int min = 0;
    for(int i=0; i<numVertices; i++)
    {
      if(S[i] == 0 && vertSet[i+1] == vertSet[k])
      {
        if(!flag)
        {
          min = f[getIndex(i+1,sigma[0])];
          flag = true;
        }
        else
        {
          if(f[getIndex(i+1,sigma[0])] < min)
          {
            min = f[getIndex(i+1,sigma[0])];
          }
        }
      }
    }

    for(int i=0; i<numVertices; i++)
    {
      if(S[i] == 0 && vertSet[i+1] == vertSet[k])
      {
        if(f[getIndex(i+1,sigma[0])] == min)
        {
          sigma[k-1] = i+1;
          int* newS = new int[numVertices];
          for(int j=0; j<numVertices; j++)
            newS[j] = S[j];
          newS[i] = 1;
          lexMinRecursion(k+1,sigma,newS,sigmaStar,f,vPart);
          delete[] newS;
        }
      }
    }
  }
  delete[] vertSet;
  return;
}


bool addFacet(Facet& f, Facet* newFacet)
{
  Facet* fPtr = &f;

  while(fPtr != NULL)
  {
    //check if *fPtr is same as newFacet
    bool same = true;
    if(newFacet->b != fPtr->b)
      same = false;
    for(int i=0; i<d && same; i++)
    {
      if(newFacet->a[i] != fPtr->a[i])
        same = false;
    }
    if(same)
      return false;
    fPtr = fPtr->next;
  }

  Facet* newF = new Facet;
  newF->a = new int[d];
  newF->next = NULL;
  for(int i=0; i<d; i++)
    newF->a[i] = newFacet->a[i];
  newF->b = newFacet->b;
  fPtr = &f;
  while(fPtr->next != NULL)
    fPtr = fPtr->next;
  fPtr->next = newF;
  return true;
}

int** getPoints(Facet f, int** points, int oldNumPoints, int& 
newNumPoints)
{
  newNumPoints = 0;
  int** newPoints = new int*[oldNumPoints];
  for(int i=0; i<oldNumPoints; i++)
  {
    int fp = 0;
    for(int j=0; j<d; j++)
    {
      fp += f.a[j]*points[i][j];
    }
/*    if(fp > f.b)
    {
      cout << "Error(2): point on wrong side of face" << endl;
      cout << "point is: ";
      for(int j=0; j<d; j++)
        cout << points[i][j] << " ";
      cout << endl << "face is: ";
      for(int j=0; j<d; j++)
        cout << f.a[j] << " ";
      cout << f.b << endl;
      exit(1);
    }*/
    if(fp == f.b)
    {
      //we have a new point
      newPoints[newNumPoints] = new int[d];
      for(int j=0; j<d; j++)
        newPoints[newNumPoints][j] = points[i][j];
      newNumPoints++;
    }
  }
  return newPoints;
}

int* rotate(int* g, int** points, int* f, int N)
{
  int xBar[d];
  int min = 0;
  int* newG = new int[d+1];

  //setup and find initial xBar
  for(int i=0; i<d; i++)
  {
    min += f[i]*points[0][i];
    xBar[i] = points[0][i];
    newG[i] = g[i];
  }
  newG[d] = g[d];
  for(int j=1; j<N; j++)
  {
    int temp = 0;
    for(int i=0; i<d; i++)
    {
      temp += f[i]*points[j][i];
    }
    if(temp < min)
    {
      min = temp;
      for(int i=0; i<d; i++)
      {
        xBar[i] = points[j][i];
      }
    }
  }


  int gxBar = 0;
  for(int i=0; i<d; i++)
  {
    gxBar += newG[i]*xBar[i];
  }

bool flag = true;
  while(gxBar != newG[d] || flag)
  {
    flag = false;
    //rotate newG
    int t1 = f[d];
    int t2 = newG[d];
    for(int i=0; i<d; i++)
    {
      t1 -= f[i]*xBar[i];
      t2 -= newG[i]*xBar[i];
    }
    for(int i=0; i<=d; i++)
    {
      newG[i] = t1*newG[i] - t2*f[i];
    }

    //update xBar
    int max = 0;
    for(int i=0; i<d; i++)
    {
      max += newG[i]*points[0][i];
      xBar[i] = points[0][i];
    }
    for(int j=1; j<N; j++)
    {
      int temp = 0;
      for(int i=0; i<d; i++)
      {
        temp += newG[i]*points[j][i];
      }
      if(temp > max)
      {
        max = temp;
        for(int i=0; i<d; i++)
        {
          xBar[i] = points[j][i];
        }
      }
    }
    gxBar = max;

  }
  return newG;
}

Facet getFacetsByPorta(int** points, int numPoints, Partition vPart, Equality* hyperP)
{
  pointListList*  pLLPtr = pLL;
  if(pLLPtr != NULL)
  {
    while(pLLPtr->next != NULL)
    {
      if(pLLPtr->numPoints == numPoints)
      {
        pointList* pL = pLLPtr->pL;
        bool flag = true;
        for(int i=0; i<numPoints && flag; i++)
        {
          for(int j=0; j<d && flag; j++)
          {
            if(points[i][j] != pL->a[j])
              flag = false;
          }
          pL = pL->next;
        }
        if(flag)
        {
          Facet initF;
          initF.a = new int[d];
          initF.next = NULL;
          Facet* p = &(pLLPtr->fac);
          while(p != NULL)
          {
            Facet newF;
            newF.a = new int[d];
            for(int i=0; i<d; i++)
              newF.a[i] = p->a[i];
            newF.b = p->b;
            newF.next = NULL;
            normalizeFacet(newF,vPart,hyperP);
            addFacet(initF,&newF);
            p = p->next;
          }
          
          return *(initF.next);
        }
      }
      pLLPtr = pLLPtr->next;
    }
  }
  Facet initF;
  initF.a = new int[d];
  initF.next = NULL;
//  Facet* fPtr = &initF;
/*  srand(time(NULL));
  int rnd = rand()%10000000;*/
  char fileName[100];
  sprintf(fileName,"tmp%d.poi",fileCount);
ofstream fout(fileName);
  if(!fout)
  {
    cout << "Error: could not open " << fileName << endl;
    exit(1);
  }
  fout << "DIM = " << d << endl;
  fout << "CONV_SECTION" << endl;
  for(int i=0; i<numPoints; i++)
  {
    for(int j=0; j<d; j++)
    {
      fout << points[i][j] << " ";
    }
    fout << endl;
  }
  fout << "END" << endl;
  fout.close();

  char cmd[100];
  sprintf(cmd,"./traf -p %s", fileName);
  time_t start,end;
  time(&start);
  system(cmd);
  time(&end);
  totSystemTime += difftime(end,start);

  //read in results of system call
  char inFileName[100];
  sprintf(inFileName,"%s.ieq",fileName);  
  ifstream fin(inFileName);
  if(!fin)
  {
    cout << "Error: could not open " << inFileName << endl;
    exit(1);
  }

    const int NN = 200;
    char s[NN];
    s[0] = 0;
    while(strchr(s,'<') == NULL)
    {
      fin.getline(s,NN);
    }
    bool flag2 = true;
    bool lastIneq = false;
    while(flag2)
    {
      if(lastIneq)
        flag2 = false;
      int* g = new int[d+1];
      for(int i=0; i<=d; i++)
        g[i] = 0;
      char* ss = strtok(s,"x");
      char* tmpS = strchr(s,')');
      ss = tmpS + sizeof(char);
      int lastCo = atoi(ss);
      if(lastCo == 0)
      {
        //ss is "-" or "+"
        bool flag = true;
        for(int i=0; i<strlen(ss) && flag; i++)
        {
          if(ss[i] == '-')
          {
            flag = false;
            lastCo = -1;
          }
          if(ss[i] == '+')
          {
            flag = false;
            lastCo = 1;
          }
        }
      }
      while(ss = strtok(NULL,"x"))
      {
        int t=0;
        while(ss[t] != '-' && ss[t] != '+' && ss[t] != 0)
          t++;
        char sss[NN];
        for(int i=0; i<t; i++)
          sss[i] = ss[i];
        sss[t] = 0;
        int ind = atoi(sss);
        g[ind-1] = lastCo;
        lastCo = atoi(ss+t*sizeof(char));
        if(lastCo == 0)
        {
          //ss is "-" or "+"
          bool flag = true;
          for(int i=t; i<strlen(ss) && flag; i++)
          {
            if(ss[i] == '-')
            {
              flag = false;
              lastCo = -1;
            }
            if(ss[i] == '+')
            {
              flag = false;
              lastCo = 1;
            }
          }
        }
      }
      int tt = 0;
      for(int i=0; i<NN; i++)
      {
        if(s[i] == '=')
        {
          tt = i+1;
          break;
        }
      }
      g[d] = atoi(s+tt*sizeof(char));

      Facet newF;
      newF.a =new int[d];
      newF.next = NULL;
      for(int i=0; i<d; i++)
        newF.a[i] = g[i];
      newF.b = g[d];
      normalizeFacet(newF,vPart,hyperP);
      addFacet(initF,&newF);
      char c;
      fin.get(c);
      if(c != '(')
      {
        flag2 = false;
      }
      else
      {
        fin.getline(s,NN);
      }
      delete[] g;
    }
    delete[] initF.a;
    if(pLLPtr != NULL)
    {
      pLLPtr->next = new pointListList;
      pLLPtr = pLLPtr->next;
    }
    else
    {
      pLL = new pointListList;
      pLLPtr = pLL;
    }
    pLLPtr->next = NULL;
    pLLPtr->numPoints = numPoints;
    pLLPtr->pL = new pointList;
    pointList* pLPtr = pLLPtr->pL;
    for(int i=0; i<numPoints; i++)
    {
      pLPtr->a = new int[d];
      for(int j=0; j<d; j++)
        pLPtr->a[j] = points[i][j];
      if(i == numPoints-1)
        pLPtr->next = NULL;
      else
      {
        pLPtr->next = new pointList;
        pLPtr = pLPtr->next;
      }
    }
//    pLLPtr->fac = *(initF.next);
    Facet* fP1 = initF.next;
    Facet* fP2 = &(pLLPtr->fac);
//    fP2 = new Facet;
    while(fP1 != NULL)
    {
      fP2->a = new int[d];
      for(int i=0; i<d; i++)
        fP2->a[i] = fP1->a[i];
      fP2->b = fP1->b;
      fP2->next = new Facet;
      fP1 = fP1->next;
      if(fP1 != NULL)
        fP2 = fP2->next;
    }
    fP2->next = NULL;
    return *(initF.next);
}

int rank(int** points, int numPoints)
{
  int** m = new int*[numPoints];
  for(int i=0; i<numPoints; i++)
  {
    m[i] = new int[d];
    for(int j=0; j<d; j++)
      m[i][j] = points[i][j];
  }
  bool foundCol = false;
  int col = -1;
  int row = -1;
  for(int i=0; i<numPoints && !foundCol; i++)
  {
    for(int j=0; j<d && !foundCol; j++)
    {
      if(m[i][j] != 0)
      {
        col = j;
        row = i;
        foundCol = true;
      }
    }
  }
  int r;
//cout << "row is: " << row << " col is: " << col << endl;
  if(col == -1 || row == -1)
    r = 0;
  //swap the first row and i
  if(row != 0)
  {
    for(int i=0; i<d; i++)
    {
      int tmp = m[row][i];
      m[row][i] = m[0][i];
      m[0][i] = tmp;
    }
  }
  for(int i=1; i<numPoints; i++)
  {
    int tmp = m[i][col];
    for(int j=col; j<d; j++)
      m[i][j] = m[i][j]*m[0][col]-m[0][j]*tmp;
  }
  r = rank2(&(m[1]),numPoints-1)+1;
  for(int i=0; i<numPoints; i++)
    delete[] m[i];
  delete[] m;
  return r;
}

int rank2(int** points, int numPoints)
{
  int** m = new int*[numPoints];
  for(int i=0; i<numPoints; i++)
  {
    m[i] = new int[d];
    for(int j=0; j<d; j++)
      m[i][j] = points[i][j];
  }
  bool foundCol = false;
  int col = -1;
  int row = -1;
  for(int j=0; j<d && !foundCol; j++)
  {
    for(int i=0; i<numPoints && !foundCol; i++)
    {
      if(m[i][j] != 0)
      {
        col = j;
        row = i;
        foundCol = true;
      }
    }
  }
//cout << "row is: " << row << " col is: " << col << endl;
  if(col == -1 || row == -1)
  {
    for(int i=0; i<numPoints; i++)
      delete[] m[i];
    delete[] m;
    return 0;
  }
  //swap the first row and i
  if(row != 0)
  {
    for(int i=0; i<d; i++)
    {
      int tmp = m[row][i];
      m[row][i] = m[0][i];
      m[0][i] = tmp;
    }
  }
  for(int i=1; i<numPoints; i++)
  {
    int tmp = m[i][col];
    for(int j=col; j<d; j++)
      m[i][j] = m[i][j]*m[0][col]-m[0][j]*tmp;
  }
  int r = rank2(&(m[1]),numPoints-1)+1;
  for(int i=0; i<numPoints; i++)
    delete[] m[i];
  delete[] m;
  return r;
}

int quotientGroupSize(int* parts)
{
  int groupSize = 1;
  for(int i=1; i<=numVertices; i++)
    groupSize = groupSize*i;

  for(int i=1; i<=numVertices; i++)
  {
    int numI = 0;
    for(int j=1; j<=numVertices; j++)
    {
      if(parts[j] == i)
        numI++;
    }
    if(numI!=0)
    {
      int tmp = 1;
      for(int j=1; j<=numI; j++)
        tmp = tmp*j;
      groupSize = groupSize/tmp;
    }
  }
  return groupSize;
}

int* getIndices(int** m, int numEq)
{
  int* numZeros = new int[numEq];
  for(int i=0; i<numEq; i++)
  {
    numZeros[i] = 0;
    for(int j=0; j<d; j++)
    {
      if(m[i][j] == 0)
        numZeros[i] = numZeros[i] + 1;
    }
  }
//  int* indices = new int[numEq];
//  for(int i=0; i<numEq; i++)
//    indices[i] = i;
  bool swapFlag = true;
  while(swapFlag)
  {
    swapFlag = false;
    for(int i=0; i<numEq-1; i++)
    {
      //check to swap i and i+1
      if(numZeros[i] < numZeros[i+1])
      {
//        int tmpI = indices[i];
//        indices[i] = indices[i+1];
//        indices[i+1] = tmpI;
        int tmpI = numZeros[i];
        numZeros[i] = numZeros[i+1];
        numZeros[i+1] = tmpI;
        for(int j=0; j<d; j++)
        {
          int tmp = m[i][j];
          m[i][j] = m[i+1][j];
          m[i+1][j] = tmp;
        }
        swapFlag = true;
      }
    }
  }
  int* rInd = new int[numEq];
  int* alreadyZeroed = new int[d];
  for(int i=0; i<numEq; i++)
    rInd[i] = -1;
  for(int i=0; i<d; i++)
    alreadyZeroed[i] = 0;
  for(int i=0; i<numEq; i++)
  {
    bool found = false;
    for(int j=0; j<d && !found; j++)
    {
      if(m[i][j] != 0 && alreadyZeroed[j] == 0)
      {
        found = true;
        rInd[i] = j;
        alreadyZeroed[j] = 1;
        for(int s=i+1; s<numEq; s++)
        {
          if(m[s][j] != 0)
          {
            int tmp1, tmp2;
            getMultipliers(m[s][j],m[i][j],tmp1,tmp2);
            for(int t=0; t<d; t++)
            {
              m[s][t] = m[s][t]*tmp2-tmp1*m[i][t];
            }
          }
        }
      }
    }
    if(!found)
    {
      cout << "Error: did not find appropriate index" << endl;
      exit(1);
    }
  }
  delete[] alreadyZeroed;
  return rInd;
}

void getMultipliers(int a, int b, int& tmp1, int& tmp2)
{
/*  if(a == 0 || b == 0)
  {
    tmp1 = 1;
    tmp2 = 1;
    return;
  }*/
  //find gcd of a and b
  int t;
  int origA = a;
  int origB = b;
  while(b != 0)
  {
    t = b;
    b = a % b;
    a = t;
  }
  //find lcm of a and b
  t = origA*origB/a;

  tmp1 = abs(t / origB);
  tmp2 = abs(t / origA);

  if(origA < 0)
  {
    tmp1 = -1*tmp1;
  }
  if(origB < 0)
  {
    tmp2 = -1*tmp2;
  }

  return;
}
