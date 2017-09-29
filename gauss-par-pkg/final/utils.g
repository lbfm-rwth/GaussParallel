BitstringToCharFct := function( P )
 local i,Chi,l,insertCT;
 insertCT := 1;
 Chi:=[];
 l := Length(P);
 for i in [1..l] do
  if P[i] = 0 then continue; fi; 
  Chi[insertCT] := i;
  insertCT := insertCT + 1;
 od;
 return Chi;
end;

CharFctToBitstring := function( l,P )
 local i,Chi;
 Chi := [];
 for i in [ 1 .. l ] do
     Chi[i] := 0;
 od;
 for i in P do 
   Chi[i]:=1;
 od;
 return Chi;
end;

REX := function( f,TT,C )
    local i,U,D,Uct,Dct,numr,row;
 
    if IsEmpty(C) then return [[],[]]; fi;
    numr := DimensionsMat(C)[1];
    Uct := 1; Dct := 1;
    U:=[]; D:=[];
    for i in [ 1 .. numr ] do
        row := C[i];
        if  i in TT then
            U[Uct] := row; Uct := Uct + 1;
        else
            D[Dct] := row; Dct := Dct + 1;
        fi;
    od;
    ConvertToMatrixRepNC( U,f );    
    ConvertToMatrixRepNC( D,f );
    return [U,D];
end;

CEX := function( f,P,C )
    local CT;
    CT := TransposedMat(C);
    CT := REX(f,P,CT);
    CT[1] := TransposedMat(CT[1]);
    CT[2] := TransposedMat(CT[2]);
    return CT;
end;

PVC := function( s,t )
 local stOrdered,u,i,l;
 stOrdered := Concatenation(s,t); 
 l := Length(s)+Length(t);
 u := [];
 for i in [ 1 .. l ] do
     u[i] := 0;
 od;
 Sort(stOrdered);
 for i in t do
  u[Position(stOrdered,i)]:=1;
 od;
  
 return [stOrdered,u];
end;

RRF := function( f,R,RR,u )
 local l,ind,indR,indRR,Rnew;
 indR:=1; indRR:=1;
 Rnew := [];
 l := Length(u);
 if R = [] then return RR; fi; if RR = [] then return R; fi;

 while (indR <= DimensionsMat(R)[1]) or (indRR <= DimensionsMat(RR)[1]) do
  ind := indR + indRR -1;
  if l < ind then
      break;
  fi;
  if u[ind] = 0 then
   Rnew[ind] := R[indR];
   indR := indR + 1;
  else
   Rnew[ind] := RR[indRR];
   indRR := indRR + 1;
  fi;
 od;
 
 ConvertToMatrixRepNC( Rnew,f );
 return Rnew;
end;

ECH := function( f,H )
 local sct,Mct,Kct,tct,Rct,EMT,m,k,M,K,R,S,N,r,s,t,i,ind,Id,one,zero,dims,dimId;

 if H = [] then return [ [],[],[],[],[] ]; fi;
 if Rank(H) = 0 then return [ [],[],[],[],[] ]; fi;   #--noetig?
 
 EMT := EchelonMatTransformation( H );
 m := TransposedMat(EMT.coeffs);
 r := TransposedMat(EMT.vectors);
 k := TransposedMat(EMT.relations);
 s := [];
 t := [];
 R := [];
 M := [];
 K := [];
 one := One(f);
 zero := Zero(f);
 Id := IdentityMat( Length(EMT.vectors),f );
 sct := 1;
 Mct := 1;
 Kct := 1;
 Rct := 1;
 tct := 1;

 ind := 1;
 dims := DimensionsMat(m);
 for i in [1..dims[1]] do
  if m[i] = zero*[1..dims[2]] then 
   s[sct] := 0;
   sct := sct + 1;
  else
   M[Mct] := m[i];
   Mct := Mct + 1;
   if not k = [] then
    K[Kct] := k[i];
    Kct := Kct + 1;
   fi;
   s[sct] := 1;
   sct := sct + 1;
  fi;
 od;
 
 ind := 1;
 dimId := DimensionsMat(Id);
 for i in [1..DimensionsMat(r)[1]] do
  if ind > dimId[1] then
   t[tct] := 0;
   tct := tct + 1;
   R[Rct] := r[i];
   Rct := Rct + 1;
   continue;
  fi;
  if r[i] = Id[ind] then
   t[tct] := 1;
   tct := tct + 1;
   ind := ind + 1;
  else
   t[tct] := 0;
   tct := tct + 1;
   R[Rct] := r[i];
   Rct := Rct + 1;
  fi;
 od; 

 M := -TransposedMat(M);
 ConvertToMatrixRepNC( M,f );
 K := TransposedMat(K);
 ConvertToMatrixRepNC( K,f );
 R := -TransposedMat(R);
 ConvertToMatrixRepNC( R,f );
 return [M,K,R,s,t];
end;
