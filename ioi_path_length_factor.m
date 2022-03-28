function pathlength = ioi_path_length_factor(lambda1, lambda2, npoints, cHBT, whichCurve)
% Compute the pathlength factor in OD x mm(-1)

% Values taken from: 
%   Kohl, Matthias,Ute Lindauer, Georg Royl, Marc K�uhl, Lorenz Gold,
%   Arno Villringer and Ulrich Dirnagl Physical model for the spectroscopic analysis of cortical
%   intrinsic optical signals,Phys. Med. Biol. 45 (2000) 3749�3764
%   figure 8a) S=50%, Tot hb 0.02
%   Values between 400-470 are set to the same value as 480.
%   Values between 690-800 are set ot the same value as 680

%  de la figure 8a) est ecrit en mm-1, mais les vraies unit�s sont des mm
%   Lambda    differential path length
%   colonne 2 OD*mm(converti par la suite en cm)

%   lambda1 Start wavelength
%   lambda2 End wavelength
%   npoints  Number of points
%   cHBT : Total hemoglobin concentration
%   whichCurve: 'Kohl' or 'Dunn'
%   deltaL: Step by half a step to evaluate values at mid point

Ext.Kohl=   ...
   [350	0.2     0.2 % longueur d'onde, sat=50% , sat=70%
    400	0.2     0.2
    410	0.2     0.2
    420	0.2     0.2
    430	0.2     0.2
    440	0.2     0.2
    450	0.2     0.2
    460	0.2     0.2
    470	0.2     0.2
    480	0.2     0.2
    490	0.2     0.2
    500	0.2     0.2
    510	0.176	0.176
    520	0.161	0.161
    530	0.125	0.125
    540	0.11	0.11
    550	0.11	0.11
    560	0.118	0.118
    570	0.118	0.118
    580	0.118	0.12272
    590	0.2     0.216
    600	0.363	0.40656
    610	0.5     0.58
    620	0.663	0.7956
    630	0.76	0.912
    640	0.89	1.068
    650	0.94	1.128
    660	1.05	1.26
    670	1.18	1.416
    680	1.25	1.5
    690	1.25	1.5
    700	1.25	1.5
    710	1.25	1.5
    720	1.25	1.5
    730	1.25	1.5
    740	1.25	1.5
    750	1.25	1.5
    760	1.25	1.5
    770	1.25	1.5
    780	1.25	1.5
    790	1.25	1.5
    820	1.25	1.5
    ];

Ext.Dunn=[450 0.52     %extrapol� Simon
        560 0.52    %Donn�es de Dunn2005
        570 0.5     %Donn�es de Dunn2005
        580 0.5     %Donn�es de Dunn2005
        590 .93     %Donn�es de Dunn2005
        600 1.58    %Donn�es de Dunn2005
        610 2.05    %Donn�es de Dunn2005
        620 2.7     %extrapol� � l'aide de Kohl
        630 3.1     %extrapol� � l'aide de Kohl
        640 3.64    %extrapol� � l'aide de Kohl
        650 3.85    %extrapol� � l'aide de Kohl
        750 3.85];  %extrapol� Simon
    
Ext.Ma = [400 0.022148 %Data from Ma et al. (2016) added 20/05/20. Were obtained from simulation that matched more closely our lens NA. (0.12 vs 0.05 for Dunn).
        402 0.01959
        404 0.016739
        406 0.013091
        408 0.009544
        410 0.007862
        412 0.006824
        414 0.00613
        416 0.006034
        418 0.006016
        420 0.006615
        422 0.007661
        424 0.009103
        426 0.010869
        428 0.012758
        430 0.014514
        432 0.016324
        434 0.021095
        436 0.025688
        438 0.030731
        440 0.042163
        442 0.051656
        444 0.071701
        446 0.087201
        448 0.121648
        450 0.173374
        452 0.225728
        454 0.288434
        456 0.322998
        458 0.347592
        460 0.377372
        462 0.412918
        464 0.433427
        466 0.467785
        468 0.499849
        470 0.526691
        472 0.554803
        474 0.580816
        476 0.604492
        478 0.626572
        480 0.64898
        482 0.666149
        484 0.67508
        486 0.683987
        488 0.692555
        490 0.697841
        492 0.705741
        494 0.714562
        496 0.723036
        498 0.730719
        500 0.73081
        502 0.731133
        504 0.727479
        506 0.731214
        508 0.721874
        510 0.712537
        512 0.700563
        514 0.685225
        516 0.664134
        518 0.625646
        520 0.58763
        522 0.544338
        524 0.49737
        526 0.45224
        528 0.410981
        530 0.371713
        532 0.338528
        534 0.315274
        536 0.295698
        538 0.282444
        540 0.272231
        542 0.269437
        544 0.272728
        546 0.281333
        548 0.296392
        550 0.316101
        552 0.336472
        554 0.355977
        556 0.374246
        558 0.383791
        560 0.392168
        562 0.396586
        564 0.390147
        566 0.373314
        568 0.3498
        570 0.324347
        572 0.299444
        574 0.280016
        576 0.271559
        578 0.278296
        580 0.306261
        582 0.354421
        584 0.434922
        586 0.542064
        588 0.677941
        590 0.847256
        592 1.040145
        594 1.252074
        596 1.481586
        598 1.702211
        600 1.998288
        602 2.156612
        604 2.342275
        606 2.51064
        608 2.645002
        610 2.794747
        612 2.940892
        614 3.100511
        616 3.205588
        618 3.305013
        620 3.410977
        622 3.506547
        624 3.602268
        626 3.693374
        628 3.771878
        630 3.846483
        632 3.924145
        634 4.000972
        636 4.053497
        638 4.096681
        640 4.140704
        642 4.185688
        644 4.231663
        646 4.274555
        648 4.311313
        650 4.348707
        652 4.386754
        654 4.425469
        656 4.464193
        658 4.502907
        660 4.535122
        662 4.565273
        664 4.595836
        666 4.626793
        668 4.657428
        670 4.687721
        672 4.718404
        674 4.748249
        676 4.775131
        678 4.801533
        680 4.82704
        682 4.852812
        684 4.878841
        686 4.905157
        688 4.928464
        690 4.948468
        692 4.964616
        694 4.980853
        696 4.995942
        698 5.009848
        700 5.023943];

if strcmp(whichCurve,'Kohl')
    E=Ext.Kohl;
    E(:,2)=E(:,2)/10; % In centimeters
elseif strcmp(whichCurve,'Dunn')
    E=Ext.Dunn;
    E(:,2)=E(:,2)/10; % In centimeters
elseif strcmp(whichCurve,'Ma')
    E=Ext.Ma;
    E(:,2) = E(:,2)/10;
else
    disp('ioi_path_length_factor: Error, no curve specified, taking data from Dunn')
    E=Ext.Dunn;
    E(:,2)=E(:,2)/10; % In centimeters
end

% Add boundaries if start and end wavelengths are outside our reach
if E(1,1)>lambda1
    E=[ lambda1 0 ;
        E(1,1)*.9999 0;
        E(:,:) ];
end

if E(end,1)<lambda2
    E=[ E(:,:);
        E(end,1)*1.00001 0;
        lambda2 0 ];
end

xi = linspace(lambda1,lambda2,npoints);
x=E(:,1); % avant x = linspace(E(1,1),E(end,1),size(E,1)); %remplac� fev 09
pathlength = E(:,2);
pathlength = interp1(x,pathlength,xi);

% trouve L pour des concentration diff�rentes de 0.02 (� partir de la
% figure 8 de kohl). On multiplie la courbe � 0.02 par un facteur de
% correction  
if strcmp(whichCurve,'Kohl') %Correction 30/07/19 ; Dunn's curve already assumes cHbT = 100uM, cHbO = 60 uM
x=[0.01 0.02 0.04 0.06 0.08 0.10 0.14 0.18]; %concentration d'h�moglobine
y=[  .33 .21 .115 .09 .075 0.062 .05 .04 ]; % valeur Da=[OD*mm] � 480 nm
correction=interp1(x,y,cHBT)/interp1(x,y,.02);
pathlength=pathlength*correction;
end
