#ifndef ROOTFILEMANAGER_H
#define ROOTFILEMANAGER_H 1

#include <Rtypes.h>
typedef Long64_t handle;
#include <map>
#include <TString.h>

class TBranch;
class TObjArray;
class TFile;
class TTree;
class TDirectory;
class TGraph;

class RootFileManager
{
  private:
    RootFileManager();
  public:
    ~RootFileManager();
    void AddFolder(TString);
    void Cd();
    void Cd(TString);
    void AddAndCd(TString);
    static RootFileManager* GetInstance() {
      return instance ? instance : instance = new RootFileManager; };
    void OpenFile( TString, const char* = "NEW" );
    void CloseFile();
    handle NewTree(TString);
    handle NewInt(handle h, TString s){return NewInt32(h,s);};
    handle NewInt32(handle, TString);
    handle NewInt64(handle, TString);
    handle NewFloat(handle h, TString s){return NewFloat32(h,s);};
    handle NewFloat32(handle, TString);
    handle NewFloat64(handle, TString);
    handle NewArray(handle, TString);
    handle NewString(handle, TString);
    handle NewGraph(handle, TString);
    void FillInt(handle h, Int_t i){FillInt32(h,i);};
    void FillInt32(handle, Int_t);
    void FillInt64(handle, Long64_t);
    void FillFloat(handle h, Float_t f){FillFloat32(h,f);};
    void FillFloat32(handle, Float_t);
    void FillFloat64(handle, Double_t);
    void FillArray(handle, TObject*);
    void FillString(handle, TString&);
    void FillGraph(handle, TGraph*);
    void AdoptTree(TTree*);
    void Adopt(TObject*);
    void CloseEntry(handle);
    void DiscardEntry(handle);
    void ClearEntries();
    void Clear();
  private:
    static RootFileManager* instance;
    handle next_handle;
    TFile* file;
    TDirectory* folder;

    std::map<handle,TTree*> t_map;
    std::map<handle,Int_t*> i_map;
    std::map<handle,Float_t*> f_map;
    std::map<handle,Long64_t*> l_map;
    std::map<handle,Double_t*> d_map;
    std::map<handle,std::pair<handle,TObjArray*> > a_map;
    std::map<handle,TString* > s_map;
    std::map<handle,TGraph** > g_map;
};

#endif

