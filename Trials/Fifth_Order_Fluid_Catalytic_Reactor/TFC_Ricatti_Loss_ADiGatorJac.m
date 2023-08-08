% This code was generated using ADiGator version 1.4
% ©2010-2014 Matthew J. Weinstein and Anil V. Rao
% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ 
% Contact: mweinstein@ufl.edu
% Bugs/suggestions may be reported to the sourceforge forums
%                    DISCLAIMER
% ADiGator is a general-purpose software distributed under the GNU General
% Public License version 3.0. While the software is distributed with the
% hope that it will be useful, both the software and generated code are
% provided 'AS IS' with NO WARRANTIES OF ANY KIND and no merchantability
% or fitness for any purpose or application.

function loss = TFC_Ricatti_Loss_ADiGatorJac(xi,setup)
global ADiGator_TFC_Ricatti_Loss_ADiGatorJac
if isempty(ADiGator_TFC_Ricatti_Loss_ADiGatorJac); ADiGator_LoadData(); end
Gator1Data = ADiGator_TFC_Ricatti_Loss_ADiGatorJac.TFC_Ricatti_Loss_ADiGatorJac.Gator1Data;
% ADiGator Start Derivative Computations
%User Line: %-------------------------------------------------------------------------%
%User Line: % %
%User Line: %
%User Line: % Authors: Kristofer Drozd
%User Line: % Created: 11/10/2021
%User Line: %
%User Line: %
%User Line: % Description:
%User Line: %
%User Line: % This function generates the loss vector for the Riccati Differential
%User Line: % Equation for the PINN/XTFC/TFC methodology.
%User Line: %
%User Line: %
%User Line: % Inputs:
%User Line: %
%User Line: % xi - A vector of orthogonal polynomial unkowns concatenated for each
%User Line: %      element in the Riccati matrix.
%User Line: %
%User Line: % setup - A structure representing inputs.
%User Line: %
%User Line: %
%User Line: % Outputs:
%User Line: %
%User Line: % loss - Loss vector.
%User Line: %
%User Line: %
%User Line: % Version:
%User Line: %
%User Line: % 11/10/2021 - Initial completion.
%User Line: %
%User Line: % %
Q = setup.ltv2DConcMats.Q;
%User Line: Q = setup.ltv2DConcMats.Q;
Sf = setup.Sf;
%User Line: Sf = setup.Sf;
A = setup.ltv2DConcMats.A;
%User Line: A = setup.ltv2DConcMats.A;
AT = setup.ltv2DConcMats.AT;
%User Line: AT = setup.ltv2DConcMats.AT;
G = setup.ltv2DConcMats.G;
%User Line: G = setup.ltv2DConcMats.G;
n = setup.tfcVariables.n;
%User Line: n = setup.tfcVariables.n;
m = setup.tfcVariables.m;
%User Line: m = setup.tfcVariables.m;
t0 = setup.t0;
%User Line: t0 = setup.t0;
tf = setup.tf;
%User Line: tf = setup.tf;
H = setup.tfcVariables.H;
%User Line: H = setup.tfcVariables.H;
Hd = setup.tfcVariables.Hd;
%User Line: Hd = setup.tfcVariables.Hd;
nx.f = size(Sf,1);
%User Line: nx = size(Sf, 1);
cada1f1 = tf - t0;
c.f = 2/cada1f1;
%User Line: c = 2/(tf-t0);
cada1f1 = n*nx.f;
CE.f = zeros(cada1f1,nx.f);
%User Line: CE = zeros(n*nx, nx);
cada1f1 = n*nx.f;
CEd.f = zeros(cada1f1,nx.f);
%User Line: CEd = zeros(n*nx, nx);
cadaforvar1.f = 1:n;
%User Line: cadaforvar1 = 1:n;
CE.dxi = zeros(6950,1);
CEd.dxi = zeros(7500,1);
for cadaforcount1 = 1:20
    k.f = cadaforvar1.f(:,cadaforcount1);
    %User Line: k = cadaforvar1(:,cadaforcount1);
    cadaforvar2.f = 1:nx.f;
    %User Line: cadaforvar2 = 1:nx;
    cada1forindex3 = reshape(Gator1Data.Index3(:,cadaforcount1),34750,5);
    cada1forindex5 = reshape(Gator1Data.Index4(:,cadaforcount1),37500,5);
    for cadaforcount2 = 1:5
        i.f = cadaforvar2.f(:,cadaforcount2);
        %User Line: i = cadaforvar2(:,cadaforcount2);
        cadaforvar3.f = 1:nx.f;
        %User Line: cadaforvar3 = 1:nx;
        cada1forindex1 = reshape(Gator1Data.Index1(:,cadaforcount2),375,5);
        cada1forindex2 = reshape(Gator1Data.Index2(:,cadaforcount2),375,5);
        cada1forindex4 = reshape(cada1forindex3(:,cadaforcount2),6950,5);
        cada1forindex6 = reshape(cada1forindex5(:,cadaforcount2),7500,5);
        for cadaforcount3 = 1:5
            j.f = cadaforvar3.f(:,cadaforcount3);
            %User Line: j = cadaforvar3(:,cadaforcount3);
            cada1f1 = i.f - 1;
            cada1f3 = cada1f1*nx.f;
            cada1f4 = cada1f3 + j.f;
            cada1f5 = cada1f4 - 1;
            cada1f7 = cada1f5*m;
            indBeg.f = cada1f7 + 1;
            %User Line: indBeg = (((i-1)*nx+j)-1)*m+1;
            cada1f1 = i.f - 1;
            cada1f3 = cada1f1*nx.f;
            cada1f4 = cada1f3 + j.f;
            cada1f5 = cada1f4 - 1;
            cada1f7 = cada1f5*m;
            indEnd.f = cada1f7 + m;
            %User Line: indEnd = (((i-1)*nx+j)-1)*m+m;
            cada1f1 = H(k.f,:);
            cada1f2 = 20;
            cada1f3 = H(cada1f2,:);
            cada1f4 = cada1f1 - cada1f3;
            cada1tempf1 = zeros(1,15);
            cada1tempf1(1:indEnd.f-indBeg.f+1) = indBeg.f:indEnd.f;
            cada1f5 = cada1tempf1;
            cada1td1 = zeros(375,1);
            cada1td1(logical(cada1forindex1(:,cadaforcount3))) = xi.dxi(nonzeros(cada1forindex1(:,cadaforcount3)));
            cada1f6dxi = cada1td1;
            cada1f6 = xi.f(cada1f5);
            cada1f7 = 15;
            cada1td1 = sparse(Gator1Data.Index12,Gator1Data.Index13,cada1f6dxi,15,375);
            cada1td1 = cada1f4*cada1td1;
            cada1td1 = cada1td1(:);
            cada1f8dxi = full(cada1td1(Gator1Data.Index14));
            cada1f8 = cada1f4*cada1f6;
            cada1f9 = Sf(i.f,j.f);
            cada1f10dxi = cada1f8dxi;
            cada1f10 = cada1f8 + cada1f9;
            cada1f11 = k.f - 1;
            cada1f13 = cada1f11*nx.f;
            cada1f14 = cada1f13 + i.f;
            CE.dxi(logical(cada1forindex4(:,cadaforcount3))) = cada1f10dxi(nonzeros(cada1forindex4(:,cadaforcount3)));
            CE.f(cada1f14,j.f) = cada1f10;
            %User Line: CE((k-1)*nx+i, j) = (H(k,:) - H(end,:))*xi(indBeg:indEnd) + Sf(i,j);
            cada1f1 = Hd(k.f,:);
            cada1f3 = c.f*cada1f1;
            cada1tempf1 = zeros(1,15);
            cada1tempf1(1:indEnd.f-indBeg.f+1) = indBeg.f:indEnd.f;
            cada1f4 = cada1tempf1;
            cada1td1 = zeros(375,1);
            cada1td1(logical(cada1forindex2(:,cadaforcount3))) = xi.dxi(nonzeros(cada1forindex2(:,cadaforcount3)));
            cada1f5dxi = cada1td1;
            cada1f5 = xi.f(cada1f4);
            cada1f6 = 15;
            cada1td1 = sparse(Gator1Data.Index15,Gator1Data.Index16,cada1f5dxi,15,375);
            cada1td1 = cada1f3*cada1td1;
            cada1td1 = cada1td1(:);
            cada1f7dxi = full(cada1td1(Gator1Data.Index17));
            cada1f7 = cada1f3*cada1f5;
            cada1f8 = k.f - 1;
            cada1f10 = cada1f8*nx.f;
            cada1f11 = cada1f10 + i.f;
            CEd.dxi(logical(cada1forindex6(:,cadaforcount3))) = cada1f7dxi(nonzeros(cada1forindex6(:,cadaforcount3)));
            CEd.f(cada1f11,j.f) = cada1f7;
            %User Line: CEd((k-1)*nx+i,j) = c*Hd(k,:)*xi(indBeg:indEnd);
        end
    end
end
cada1f1 = n*nx.f;
cada1f2 = cada1f1*nx.f;
L.f = zeros(cada1f2,1);
%User Line: L = zeros(n*nx*nx, 1);
cnt.f =  0;
%User Line: cnt = 0;
cadaforvar4.f = 1:n;
%User Line: cadaforvar4 = 1:n;
L.dxi = zeros(63100,1);
for cadaforcount4 = 1:20
    k.f = cadaforvar4.f(:,cadaforcount4);
    %User Line: k = cadaforvar4(:,cadaforcount4);
    cadaforvar5.f = 1:nx.f;
    %User Line: cadaforvar5 = 1:nx;
    cada1forindex1 = reshape(Gator1Data.Index5(:,cadaforcount4),1875,5);
    cada1forindex3 = reshape(Gator1Data.Index6(:,cadaforcount4),375,5);
    cada1forindex5 = reshape(Gator1Data.Index7(:,cadaforcount4),375,5);
    cada1forindex6 = reshape(Gator1Data.Index8(:,cadaforcount4),375,5);
    cada1forindex8 = reshape(Gator1Data.Index9(:,cadaforcount4),1875,5);
    cada1forindex10 = reshape(Gator1Data.Index10(:,cadaforcount4),315500,5);
    for cadaforcount5 = 1:5
        i.f = cadaforvar5.f(:,cadaforcount5);
        %User Line: i = cadaforvar5(:,cadaforcount5);
        cadaforvar6.f = 1:nx.f;
        %User Line: cadaforvar6 = 1:nx;
        cada1forindex2 = reshape(cada1forindex1(:,cadaforcount5),375,5);
        cada1forindex4 = cada1forindex3(:,cadaforcount5);
        cada1forindex7 = cada1forindex6(:,cadaforcount5);
        cada1forindex11 = reshape(cada1forindex10(:,cadaforcount5),63100,5);
        for cadaforcount6 = 1:5
            j.f = cadaforvar6.f(:,cadaforcount6);
            %User Line: j = cadaforvar6(:,cadaforcount6);
            cnt.f = cnt.f + 1;
            %User Line: cnt = cnt+1;
            S.f = zeros(5,1);
            S.dxi = zeros(1875,1);
            %User Line: S = zeros(nx,1);
            cadaforvar7.f = 1:nx.f;
            %User Line: cadaforvar7 = 1:nx;
            cada1forindex9 = reshape(cada1forindex8(:,cadaforcount6),375,5);
            for cadaforcount7 = 1:5
                z.f = cadaforvar7.f(:,cadaforcount7);
                %User Line: z = cadaforvar7(:,cadaforcount7);
                cada1f1 = k.f - 1;
                cada1f3 = cada1f1*nx.f;
                cada1f4 = cada1f3 + i.f;
                cada1td1 = zeros(375,1);
                cada1td1(logical(cada1forindex7)) = CE.dxi(nonzeros(cada1forindex7));
                cada1f5dxi = cada1td1;
                cada1f5 = CE.f(cada1f4,:);
                cada1f6 = k.f - 1;
                cada1f8 = cada1f6*nx.f;
                cada1f9 = cada1f8 + 1;
                cada1f10 = k.f - 1;
                cada1f12 = cada1f10*nx.f;
                cada1f13 = cada1f12 + nx.f;
                cada1tempf1 = zeros(1,5);
                cada1tempf1(1:cada1f13-cada1f9+1) = cada1f9:cada1f13;
                cada1f14 = cada1tempf1;
                cada1f15 = G(cada1f14,z.f);
                cada1f16 = 5;
                cada1td1 = sparse(Gator1Data.Index18,Gator1Data.Index19,cada1f5dxi,5,375);
                cada1td1 = cada1f15.'*cada1td1;
                cada1td1 = cada1td1(:);
                cada1f17dxi = full(cada1td1(Gator1Data.Index20));
                cada1f17 = cada1f5*cada1f15;
                cada1f18 = k.f - 1;
                cada1f20 = cada1f18*nx.f;
                cada1f21 = cada1f20 + z.f;
                cada1td1 = zeros(375,1);
                cada1td1(logical(cada1forindex9(:,cadaforcount7))) = CE.dxi(nonzeros(cada1forindex9(:,cadaforcount7)));
                cada1f22dxi = cada1td1;
                cada1f22 = CE.f(cada1f21,j.f);
                cada1td1 = cada1f22.*cada1f17dxi;
                cada1td1 = cada1td1 + cada1f17.*cada1f22dxi;
                cada1f24dxi = cada1td1;
                cada1f24 = cada1f17*cada1f22;
                S.dxi(logical(Gator1Data.Index11(:,cadaforcount7))) = cada1f24dxi(nonzeros(Gator1Data.Index11(:,cadaforcount7)));
                S.f(z.f) = cada1f24;
                %User Line: S(z) = CE((k-1)*nx+i,:)*G((k-1)*nx+1:(k-1)*nx+nx,z)*CE((k-1)*nx+z,j);
            end
            cada1f1 = k.f - 1;
            cada1f3 = cada1f1*nx.f;
            cada1f4 = cada1f3 + i.f;
            cada1td1 = zeros(375,1);
            cada1td1(logical(cada1forindex2(:,cadaforcount6))) = CEd.dxi(nonzeros(cada1forindex2(:,cadaforcount6)));
            cada1f5dxi = cada1td1;
            cada1f5 = CEd.f(cada1f4,j.f);
            cada1f6 = k.f - 1;
            cada1f8 = cada1f6*nx.f;
            cada1f9 = cada1f8 + i.f;
            cada1td1 = zeros(375,1);
            cada1td1(logical(cada1forindex4)) = CE.dxi(nonzeros(cada1forindex4));
            cada1f10dxi = cada1td1;
            cada1f10 = CE.f(cada1f9,:);
            cada1f11 = k.f - 1;
            cada1f13 = cada1f11*nx.f;
            cada1f14 = cada1f13 + 1;
            cada1f15 = k.f - 1;
            cada1f17 = cada1f15*nx.f;
            cada1f18 = cada1f17 + nx.f;
            cada1tempf1 = zeros(1,5);
            cada1tempf1(1:cada1f18-cada1f14+1) = cada1f14:cada1f18;
            cada1f19 = cada1tempf1;
            cada1f20 = A(cada1f19,j.f);
            cada1f21 = 5;
            cada1td1 = sparse(Gator1Data.Index21,Gator1Data.Index22,cada1f10dxi,5,375);
            cada1td1 = cada1f20.'*cada1td1;
            cada1td1 = cada1td1(:);
            cada1f22dxi = full(cada1td1(Gator1Data.Index23));
            cada1f22 = cada1f10*cada1f20;
            cada1td1 = cada1f5dxi;
            cada1td1 = cada1td1 + cada1f22dxi;
            cada1f23dxi = cada1td1;
            cada1f23 = cada1f5 + cada1f22;
            cada1f24 = k.f - 1;
            cada1f26 = cada1f24*nx.f;
            cada1f27 = cada1f26 + i.f;
            cada1f28 = AT(cada1f27,:);
            cada1f29 = k.f - 1;
            cada1f31 = cada1f29*nx.f;
            cada1f32 = cada1f31 + 1;
            cada1f33 = k.f - 1;
            cada1f35 = cada1f33*nx.f;
            cada1f36 = cada1f35 + nx.f;
            cada1tempf1 = zeros(1,5);
            cada1tempf1(1:cada1f36-cada1f32+1) = cada1f32:cada1f36;
            cada1f37 = cada1tempf1;
            cada1td1 = zeros(375,1);
            cada1td1(logical(cada1forindex5(:,cadaforcount6))) = CE.dxi(nonzeros(cada1forindex5(:,cadaforcount6)));
            cada1f38dxi = cada1td1;
            cada1f38 = CE.f(cada1f37,j.f);
            cada1f39 = 5;
            cada1td1 = sparse(Gator1Data.Index24,Gator1Data.Index25,cada1f38dxi,5,375);
            cada1td1 = cada1f28*cada1td1;
            cada1td1 = cada1td1(:);
            cada1f40dxi = full(cada1td1(Gator1Data.Index26));
            cada1f40 = cada1f28*cada1f38;
            cada1td1 = cada1f23dxi;
            cada1td1 = cada1td1 + cada1f40dxi;
            cada1f41dxi = cada1td1;
            cada1f41 = cada1f23 + cada1f40;
            cada1f42 = 5;
            cada1td1 = zeros(5,375);
            cada1td1(Gator1Data.Index27) = S.dxi;
            cada1td1 = sum(cada1td1,1);
            cada1f43dxi = cada1td1(:);
            cada1f43 = sum(S.f);
            cada1td1 = cada1f41dxi;
            cada1td1 = cada1td1 + -cada1f43dxi;
            cada1f44dxi = cada1td1;
            cada1f44 = cada1f41 - cada1f43;
            cada1f45 = k.f - 1;
            cada1f47 = cada1f45*nx.f;
            cada1f48 = cada1f47 + i.f;
            cada1f49 = Q(cada1f48,j.f);
            cada1f50dxi = cada1f44dxi;
            cada1f50 = cada1f44 + cada1f49;
            L.dxi(logical(cada1forindex11(:,cadaforcount6))) = cada1f50dxi(nonzeros(cada1forindex11(:,cadaforcount6)));
            L.f(cnt.f) = cada1f50;
            %User Line: L(cnt) = CEd((k-1)*nx+i,j) + CE((k-1)*nx+i,:)*A((k-1)*nx+1:(k-1)*nx+nx,j) + AT((k-1)*nx+i,:)*CE((k-1)*nx+1:(k-1)*nx+nx,j) - sum(S) + Q((k-1)*nx+i,j);
        end
    end
end
loss.dxi = L.dxi; loss.f = L.f;
%User Line: loss = L;
loss.dxi_size = [500,375];
loss.dxi_location = Gator1Data.Index28;
end


function ADiGator_LoadData()
global ADiGator_TFC_Ricatti_Loss_ADiGatorJac
ADiGator_TFC_Ricatti_Loss_ADiGatorJac = load('TFC_Ricatti_Loss_ADiGatorJac.mat');
return
end