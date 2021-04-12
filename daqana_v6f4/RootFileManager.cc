#include "RootFileManager.hh"

#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TObjArray.h>
#include <TTree.h>

#include <stdexcept>

RootFileManager* RootFileManager::instance(0);

RootFileManager::RootFileManager() : next_handle(0), file(0), folder(0)
{
  if( instance )
    throw std::runtime_error("Trying to create more than one instance of RootFileManager.");
  instance = this;
}

RootFileManager::~RootFileManager()
{
  Clear();
}

void RootFileManager::AddFolder(TString fname)
{
  if( !file || !file->IsOpen() )
    throw std::runtime_error("File not open.");
  folder->mkdir(fname);
}

void RootFileManager::Cd()
{
  if( !file || !file->IsOpen() )
    throw std::runtime_error("File not open.");
  file->cd();
  folder = file->CurrentDirectory();
}

void RootFileManager::Cd(TString fname)
{
  if( !file || !file->IsOpen() )
    throw std::runtime_error("File not open.");
  file->Cd(fname);
  folder = file->CurrentDirectory();
}

void RootFileManager::AddAndCd(TString fname)
{
  if( !file || !file->IsOpen() )
    throw std::runtime_error("File not open.");
  if( !folder->GetDirectory(fname) )
    AddFolder(fname);
  file->Cd(fname);
  folder = file->CurrentDirectory();
}

void RootFileManager::OpenFile(TString filename, const char* opt)
{
  if( file )
    throw std::runtime_error("Current file not closed when trying to open new file.");
  file = new TFile(filename, opt);
  folder = file;
  if( !file->IsOpen() )
    throw std::runtime_error("Could not open root file for writing.");
}

void RootFileManager::CloseFile()
{
  if( !file || !file->IsOpen() )
    throw std::runtime_error("File already closed.");
  file->cd();

  std::map<handle, TTree*>::iterator t_itr = t_map.begin();
  while( t_itr != t_map.end() )
    (t_itr++)->second->Write();

  file->Close();
  delete file;
  folder = 0;
  file = 0;
}

handle RootFileManager::NewTree(TString name)
{
  t_map[next_handle] = new TTree(name, name);
  return next_handle++;
}

handle RootFileManager::NewInt32(handle hndl, TString name)
{
  Int_t* ptr = new Int_t;
  *ptr = 0;
  t_map[hndl]->Branch(name,ptr,name+"/I",64000);
  i_map[next_handle] = ptr;
  return next_handle++;
}

handle RootFileManager::NewFloat32(handle hndl, TString name)
{
  Float_t* ptr = new Float_t;
  *ptr = 0;
  t_map[hndl]->Branch(name,ptr,name+"/F",64000);
  f_map[next_handle] = ptr;
  return next_handle++;
}

handle RootFileManager::NewInt64(handle hndl, TString name)
{
  Long64_t* ptr = new Long64_t;
  *ptr = 0;
  t_map[hndl]->Branch(name,ptr,name+"/L",64000);
  l_map[next_handle] = ptr;
  return next_handle++;
}

handle RootFileManager::NewFloat64(handle hndl, TString name)
{
  Double_t* ptr = new Double_t;
  *ptr = 0;
  t_map[hndl]->Branch(name,ptr,name+"/D",64000);
  d_map[next_handle] = ptr;
  return next_handle++;
}

handle RootFileManager::NewArray(handle hndl, TString name)
{
  TObjArray* ptr = new TObjArray;
  ptr->SetOwner(true); // Will delete objects when clearing
  t_map[hndl]->Branch(name, "TObjArray", ptr,64000);
  a_map[next_handle] = std::pair<handle, TObjArray*>(hndl,ptr);
  return next_handle++;
}

handle RootFileManager::NewString(handle hndl, TString name)
{
  TString* ptr = new TString;
  t_map[hndl]->Branch(name,"TString", ptr,64000);
  s_map[next_handle] = ptr;
  return next_handle++;
}

handle RootFileManager::NewGraph(handle hndl, TString name)
{
  TGraph** ptr = new TGraph*;
  *ptr = new TGraph;
  t_map[hndl]->Branch(name,"TGraph", ptr,64000);
  g_map[next_handle] = ptr;
  return next_handle++;
}

void RootFileManager::FillInt32(handle hndl, Int_t data)
{
  *i_map[hndl] = data;
}

void RootFileManager::FillFloat32(handle hndl, Float_t data)
{
  *f_map[hndl] = data;
}

void RootFileManager::FillInt64(handle hndl, Long64_t data)
{
  *l_map[hndl] = data;
}

void RootFileManager::FillFloat64(handle hndl, Double_t data)
{
  *d_map[hndl] = data;
}

void RootFileManager::FillArray(handle hndl, TObject* data)
{
  a_map[hndl].second->Add(data);
}

void RootFileManager::FillString(handle hndl, TString& data)
{
  *s_map[hndl] = data;
}

void RootFileManager::FillGraph(handle hndl, TGraph* data)
{
  TGraph* ptr = data;
  TGraph** pptr = g_map[hndl];
  delete *pptr;
  *pptr = ptr;
}

void RootFileManager::AdoptTree( TTree* t )
{
  folder->cd();
  t->CloneTree()->Write();
}

void RootFileManager::Adopt( TObject* o )
{
  folder->cd();
  o->Write();
  delete o;
}

void RootFileManager::CloseEntry(handle hndl)
{
  t_map[hndl]->Fill();
  std::map<handle, std::pair<handle,TObjArray*> >::iterator itr = a_map.begin();
  while( itr != a_map.end() )
  {
    if( itr->second.first == hndl )
      itr->second.second->Delete();
    ++itr;
  }
}

void RootFileManager::ClearEntries()
{
  for( std::map<handle, std::pair<handle,TObjArray*> >::iterator itr = a_map.begin();
      itr != a_map.end(); ++itr )
    itr->second.second->Delete();
  for( std::map<handle, Int_t* >::iterator itr = i_map.begin();
      itr != i_map.end(); ++itr ) *itr->second = 0;
  for( std::map<handle, Float_t* >::iterator itr = f_map.begin();
      itr != f_map.end(); ++itr ) *itr->second = 0;
  for( std::map<handle, Long64_t* >::iterator itr = l_map.begin();
      itr != l_map.end(); ++itr ) *itr->second = 0;
  for( std::map<handle, Double_t* >::iterator itr = d_map.begin();
      itr != d_map.end(); ++itr ) *itr->second = 0;
  for( std::map<handle, TString* >::iterator itr = s_map.begin();
      itr != s_map.end(); ++itr ) *itr->second = "";
}

void RootFileManager::DiscardEntry(handle hndl)
{
  for( std::map<handle, std::pair<handle,TObjArray*> >::iterator itr = a_map.begin();
      itr != a_map.end(); ++itr )
    if( itr->second.first == hndl )
      itr->second.second->Delete();
}

// Are we leaking TTrees and adopted objects? Probably...
void RootFileManager::Clear()
{
  if( file && file->IsOpen() ) CloseFile();
  next_handle = 0;
  t_map.clear();
  i_map.clear();
  l_map.clear();
  f_map.clear();
  d_map.clear();
  std::map<handle, std::pair<handle,TObjArray*> >::iterator itr = a_map.begin();
  while( itr != a_map.end() )
  {
    TObjArray* ptr = (itr++)->second.second;
    ptr->Delete();
    delete ptr;
  }
  a_map.clear();
  folder = 0;
}

