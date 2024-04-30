within ;
package TJUthermo
  package Interface
    partial connector Flange
    Modelica.SIunits.Pressure p;
    flow Modelica.SIunits.MassFlowRate m_flow;
    stream Modelica.SIunits.SpecificEnthalpy h;//[j/kg]
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Flange;

   connector Flange_a
    extends Flange;
     annotation (
       Icon(graphics={Ellipse(
       extent={{-100,100},{100,-100}},
       lineColor={0,0,255},
       fillColor={255,0,0},
       fillPattern=FillPattern.Solid)}), uses(Modelica(version="3.2.1")));
   end Flange_a;

  connector Flange_b
    extends Flange;
     annotation (
       Icon(graphics={Ellipse(
       extent={{-100,100},{100,-100}},
       lineColor={0,0,255},
       fillColor={255,255,0},
       fillPattern=FillPattern.Solid)}), uses(Modelica(version="3.2.1")));
  end Flange_b;

  connector ThermoT
  Modelica.SIunits.Temperature T;
  flow Modelica.SIunits.HeatFlowRate Q;
      annotation (Icon(graphics={Rectangle(
              extent={{-68,24},{74,-2}},
              lineColor={28,108,200},
              fillColor={95,95,95},
              fillPattern=FillPattern.Solid), Text(
              extent={{-46,-10},{46,-38}},
              lineColor={28,108,200},
              fillColor={95,95,95},
              fillPattern=FillPattern.Solid,
              textString="wall temperature")}));
  end ThermoT;

  model single_wall_port

    ThermoT tin
      annotation (Placement(transformation(extent={{-10,-52},{10,-32}})));
    ThermoT tout
      annotation (Placement(transformation(extent={{-12,28},{8,48}})));
     parameter Modelica.SIunits.Mass M=0.1;
     parameter Modelica.SIunits.SpecificHeatCapacity Cp;
     Modelica.SIunits.Temperature T(start=313);
  equation
    M*Cp*der(T)=tout.Q+tin.Q;
    T=tin.T;
    T=tout.T;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
            extent={{-70,36},{78,-36}},
            lineColor={28,108,200},
            fillColor={255,170,213},
            fillPattern=FillPattern.Solid), Text(
            extent={{-52,26},{58,-22}},
            lineColor={28,108,200},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            textString="Wall")}), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end single_wall_port;

  connector mutiThermoP
  parameter Integer N(min=1)=2 "Number of nodes";
    Modelica.SIunits.Temperature T[N] "Temperature at the nodes";
    flow Modelica.SIunits.HeatFlowRate Q[N] "Heat flux at the nodes";
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
            extent={{-76,22},{80,12}},
            lineColor={28,108,200},
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-76,-14},{80,-26}},
            lineColor={28,108,200},
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{10,12},{-10,-14}},
            lineColor={28,108,200},
            fillColor={170,255,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{10,10},{66,-12}},
            lineColor={28,108,200},
            fillColor={170,255,255},
            fillPattern=FillPattern.Solid,
            textString="-wall"),
          Text(
            extent={{-72,10},{-20,-12}},
            lineColor={28,108,200},
            fillColor={170,255,255},
            fillPattern=FillPattern.Solid,
            textString="Thermo")}),                                Diagram(
          coordinateSystem(preserveAspectRatio=false)));

  end mutiThermoP;

  model ThermoConvet
    parameter Integer N=10;
    ThermoT single[N];
    mutiThermoP multi(final N=N);
  equation
    for i in 1:N loop
      single[i].T=multi.T[i];
      single[i].Q=-multi.Q[i];
    end for;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end ThermoConvet;

    partial connector SimpleFlange
    flow Modelica.SIunits.MassFlowRate m_flow;
     Modelica.SIunits.Temperature T; //[j/kg]
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));

    end SimpleFlange;

    connector SimpleFlange_a
      extends SimpleFlange;
      annotation (Icon(graphics={Ellipse(
              extent={{-98,86},{98,-86}},
              lineColor={28,108,200},
              fillColor={255,0,0},
              fillPattern=FillPattern.Forward)}));
    end SimpleFlange_a;

    connector SimpleFlange_b
    extends Interface.SimpleFlange;
    annotation (Icon(graphics={Ellipse(
            extent={{100,82},{-100,-94}},
            lineColor={28,108,200},
            fillColor={255,255,85},
            fillPattern=FillPattern.Forward)}));
    end SimpleFlange_b;

    model DeMix_node
     parameter Integer N=3;
      SimpleFlange_a multi[N];
      Modelica.SIunits.Temperature T;
      Modelica.SIunits.MassFlowRate m_flow;
      SimpleFlange_a single
        annotation (Placement(transformation(extent={{-66,-8},{-46,12}})));
    equation
      m_flow=single.m_flow/N;
      T=single.T;
      for i in 1:N loop
       single.T=multi[i].T;
       m_flow=-multi[i].m_flow;
      end for
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end DeMix_node;

    connector mutiNode_liquid
    parameter Integer N(min=1)=2 "Number of nodes";
      Modelica.SIunits.Temperature T[N] "Temperature at the nodes";
      flow Modelica.SIunits.MassFlowRate m_flow[N] "Heat flux at the nodes";
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end mutiNode_liquid;

    connector SimpleFlangeMulti
      parameter Integer N(min=1)=2;
    flow Modelica.SIunits.MassFlowRate m_flow[N];
     Modelica.SIunits.Temperature T[N]; //[j/kg]
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Ellipse(
              extent={{-94,98},{98,-96}},
              lineColor={28,108,200},
              fillColor={175,175,175},
              fillPattern=FillPattern.Solid), Text(
              extent={{-48,42},{48,-26}},
              lineColor={28,108,200},
              fillColor={175,175,175},
              fillPattern=FillPattern.Solid,
              textString="M")}),                                     Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end SimpleFlangeMulti;

    model nodeConvet_liquid
     parameter Integer N=10;
      SimpleFlange_a single[N];
      SimpleFlangeMulti multi(final N=N);
    equation
      for i in 1:N loop
        single[i].T=multi.T[i];
        single[i].m_flow=-multi.m_flow[i];
      end for;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end nodeConvet_liquid;
  end Interface;

  partial package Media
     extends ExternalMedia.Media.CoolPropMedium;

  package R245fa_CP
      extends ExternalMedia.Media.CoolPropMedium(mediumName="R245fa");
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
  end R245fa_CP;

    package R123_CP
      extends ExternalMedia.Media.CoolPropMedium(mediumName="R245fa");
    end R123_CP;

    package MediaUSD

      model R125_USD
        //parameter Modelica.SIunits.Temperature T=300;
        //parameter Modelica.SIunits.Density rou=89;
        constant Modelica.SIunits.Temperature Tc=339.173 "临界温度";
        constant Modelica.SIunits.Density rouc=573.58 "临界密度";
        constant Modelica.SIunits.Pressure Pc=3.6177e6 "临界压力";
        constant Modelica.SIunits.MolarMass MM=120.02e-3 "摩尔质量";
        constant TJUthermo.DataRecord.GasConstant R=8.314472
          "气体摩尔常数";
        //--------------------系数-----------------------//
        model R125state_ph
        Modelica.SIunits.Temperature T(start=300);
        Modelica.SIunits.Density  rou(start=1000);
          parameter Modelica.SIunits.SpecificEnthalpy h;
        parameter Modelica.SIunits.Pressure P;
        protected
        Real A;
        Real B;
        Real C;
        Real tao;
        Real delta;
        //Real x;
        //TJUthermo.DataRecord.SatProperty Sat;
        parameter Real a[6]={37.2674,8.88404,-49.8651,2.303,5.086,7.3};
        parameter Real b[6]={0,0,0,0.92578,2.22895,5.03283};
        parameter Real N[18]={5.280760000,-8.676580000,0.750112700,0.759002300,0.014518990,4.777189000,-3.330988000,3.775673000,
        -2.290919000,0.888826800,-0.623486400,-0.041272630,-0.084553890,-0.130875200,0.008344962,-1.532005000,-0.058836490,0.022966580};
        parameter Real d[18]={1,1,1,2,4,1,1,2,2,3,4,5,1,5,1,2,3,5};
        parameter Real t[18]={0.669,1.050,2.750,0.956,1.000,2.000,2.750,2.380,3.370,3.470,2.630,3.450,0.720,4.230,0.200,4.500,29.000,24.000};
        parameter Real l[18]={0,0,0,0,0,1,1,1,1,1,1,1,2,2,3,2,3,3};
        parameter Real m[18]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.7,7,6};
        equation
         // if T<Tc then
          // Sat=TJUthermo.Media.MediaUSD.R125_USD.Sat_property_T(T);
          //end if;
        tao=Tc/T;
        delta=rou/rouc;
        algorithm
        A:=0;
        for i in 1:5 loop
            A := A + d[i]*N[i]*delta^d[i]*tao^t[i];
        end for;
        for i in 6:15 loop
            A := A + N[i]*delta^d[i]*tao^t[i]*exp(-delta^l[i])*(d[i] - l[i]*delta^l[i]);
          end for;
        for i in 16:18 loop
            A := A + N[i]*delta^d[i]*tao^t[i]*exp(-delta^l[i])*exp(-tao^m[i])*(d[i] -
              l[i]*delta^l[i]);
          end for;
          C := 0;
        for i in 1:5 loop
         C:=C+t[i]*N[i]*delta^d[i]*tao^t[i];
          end for;
        for i in 6:15 loop
         C:=C+t[i]*N[i]*delta^d[i]*tao^t[i]*exp(-delta^l[i]);
          end for;
        for  i in 16:18 loop
            C := C + N[i]*delta^d[i]*tao^t[i]*exp(-delta^l[i])*exp(-tao^m[i])*(t[i] -
              m[i]*tao^m[i]);
          end for;
        //P:=rou/(MM)*R*T*(1+A);
        B:=(a[2]+(-0.1)*tao^(-0.1)+a[4]*b[4]/(exp(b[4]*tao)-1)+a[5]*b[5]/(exp(b[5]*tao)-1)+a[6]*b[6]/(exp(b[6]*tao)-1))*tao+1;
        equation
        //P=rou/(MM)*R*T*(1+A);
          // if T<Tc and rou<Sat.roul and rou>Sat.rouv then
          //  x=(Sat.rouv*Sat.roul/rou-Sat.rouv)/(Sat.roul-Sat.rouv);
          //  h=x*Sat.hv+(1-x)*Sat.hl;
          //  P=Sat.Psat;
         // else
            h=R*T*(B+C+A+1)/MM+6.3121e4;
            P=rou/(MM)*R*T*(1+A);
            //x:=2;
          //  end if;
        end R125state_ph;

        function Sat_property_T
              input Modelica.SIunits.Temperature T;
              output TJUthermo.DataRecord.SatProperty Sat;
        protected
        parameter Real Nl1=1.6684;
        parameter Real Nl2=0.88415;
        parameter Real Nl3=0.44383;
        parameter Real Ns1=-2.8403;
        parameter Real Ns2=-7.2738;
        parameter Real Ns3=-21.890;
        parameter Real Ns4=-58.825;
        parameter Real Np1=-7.5295;
        parameter Real Np2=1.9026;
        parameter Real Np3=-2.2966;
        parameter Real Np4=-3.4480;
        Real theta;
        algorithm
        //---------------------求饱和点物性密度------------------------//
        theta:=1-T/Tc;
        Sat.roul:=(1+Nl1*theta^(1/3)+Nl2*theta^0.6+Nl3*theta^2.9)*rouc;
        Sat.rouv:=exp(Ns1*theta^0.38+Ns2*theta^1.22+Ns3*theta^3.3+Ns4*theta^6.9)*rouc;
        //------------------------饱和温度下压力-----------------------------//
        Sat.Psat:=Pc*exp(Tc/T*(Np1*theta + Np2*theta^1.5 + Np3*theta^2.3 + Np4*theta^4.6));
        //-----饱和温度下焓值----//
        Sat.hl:=TJUthermo.Media.MediaUSD.R125_USD.SimpleEnthalpy_Trou(T,Sat.roul);
        Sat.hv:=TJUthermo.Media.MediaUSD.R125_USD.SimpleEnthalpy_Trou(T,Sat.rouv);
        end Sat_property_T;
        //-----根据温度密度求焓值----//
        function SpecificEnthalpy_Trou
          input Modelica.SIunits.Temperature T;
          input Modelica.SIunits.Density  rou;
          output Modelica.SIunits.SpecificEnthalpy h;
        protected
        Real A;
        Real B;
        Real C;
        Real tao;
        Real delta;
        Real x;
        TJUthermo.DataRecord.SatProperty Sat;
        parameter Real a[6]={37.2674,8.88404,-49.8651,2.303,5.086,7.3};
        parameter Real b[6]={0,0,0,0.92578,2.22895,5.03283};
        parameter Real N[18]={5.280760000,-8.676580000,0.750112700,0.759002300,0.014518990,4.777189000,-3.330988000,3.775673000,
        -2.290919000,0.888826800,-0.623486400,-0.041272630,-0.084553890,-0.130875200,0.008344962,-1.532005000,-0.058836490,0.022966580};
        parameter Real d[18]={1,1,1,2,4,1,1,2,2,3,4,5,1,5,1,2,3,5};
        parameter Real t[18]={0.669,1.050,2.750,0.956,1.000,2.000,2.750,2.380,3.370,3.470,2.630,3.450,0.720,4.230,0.200,4.500,29.000,24.000};
        parameter Real l[18]={0,0,0,0,0,1,1,1,1,1,1,1,2,2,3,2,3,3};
        parameter Real m[18]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.7,7,6};
        algorithm
          if T<Tc then
           Sat:=TJUthermo.Media.MediaUSD.R125_USD.Sat_property_T(T);
          end if;
        tao:=Tc/T;
        delta:=rou/rouc;
        A:=0;
        for i in 1:5 loop
            A := A + d[i]*N[i]*delta^d[i]*tao^t[i];
        end for;
        for i in 6:15 loop
            A := A + N[i]*delta^d[i]*tao^t[i]*exp(-delta^l[i])*(d[i] - l[i]*delta^l[i]);
          end for;
        for i in 16:18 loop
            A := A + N[i]*delta^d[i]*tao^t[i]*exp(-delta^l[i])*exp(-tao^m[i])*(d[i] -
              l[i]*delta^l[i]);
          end for;
          C := 0;
        for i in 1:5 loop
         C:=C+t[i]*N[i]*delta^d[i]*tao^t[i];
          end for;
        for i in 6:15 loop
         C:=C+t[i]*N[i]*delta^d[i]*tao^t[i]*exp(-delta^l[i]);
          end for;
        for  i in 16:18 loop
            C := C + N[i]*delta^d[i]*tao^t[i]*exp(-delta^l[i])*exp(-tao^m[i])*(t[i] -
              m[i]*tao^m[i]);
          end for;
        //P:=rou/(MM)*R*T*(1+A);
        B:=(a[2]+(-0.1)*tao^(-0.1)+a[4]*b[4]/(exp(b[4]*tao)-1)+a[5]*b[5]/(exp(b[5]*tao)-1)+a[6]*b[6]/(exp(b[6]*tao)-1))*tao+1;
        if T<Tc and rou<Sat.roul and rou>Sat.rouv then
            x:=(Sat.rouv*Sat.roul/rou-Sat.rouv)/(Sat.roul-Sat.rouv);
            h:=x*Sat.hv+(1-x)*Sat.hl;
          else
            h:=R*T*(B+C+A+1)/MM+6.3121e4;
            //x:=2;
        end if;
        end SpecificEnthalpy_Trou;

        function SimpleEnthalpy_Trou
             input Modelica.SIunits.Temperature T;
          input Modelica.SIunits.Density  rou;
          output Modelica.SIunits.SpecificEnthalpy h;
        protected
        Real A;
        Real B;
        Real C;
        Real tao;
        Real delta;
        parameter Real a[6]={37.2674,8.88404,-49.8651,2.303,5.086,7.3};
        parameter Real b[6]={0,0,0,0.92578,2.22895,5.03283};
        parameter Real N[18]={5.280760000,-8.676580000,0.750112700,0.759002300,0.014518990,4.777189000,-3.330988000,3.775673000,
        -2.290919000,0.888826800,-0.623486400,-0.041272630,-0.084553890,-0.130875200,0.008344962,-1.532005000,-0.058836490,0.022966580};
        parameter Real d[18]={1,1,1,2,4,1,1,2,2,3,4,5,1,5,1,2,3,5};
        parameter Real t[18]={0.669,1.050,2.750,0.956,1.000,2.000,2.750,2.380,3.370,3.470,2.630,3.450,0.720,4.230,0.200,4.500,29.000,24.000};
        parameter Real l[18]={0,0,0,0,0,1,1,1,1,1,1,1,2,2,3,2,3,3};
        parameter Real m[18]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.7,7,6};
        algorithm
        tao:=Tc/T;
        delta:=rou/rouc;
        A:=0;
        for i in 1:5 loop
            A := A + d[i]*N[i]*delta^d[i]*tao^t[i];
        end for;
        for i in 6:15 loop
            A := A + N[i]*delta^d[i]*tao^t[i]*exp(-delta^l[i])*(d[i] - l[i]*delta^l[i]);
          end for;
        for i in 16:18 loop
            A := A + N[i]*delta^d[i]*tao^t[i]*exp(-delta^l[i])*exp(-tao^m[i])*(d[i] -
              l[i]*delta^l[i]);
          end for;
          C := 0;
        for i in 1:5 loop
         C:=C+t[i]*N[i]*delta^d[i]*tao^t[i];
          end for;
        for i in 6:15 loop
         C:=C+t[i]*N[i]*delta^d[i]*tao^t[i]*exp(-delta^l[i]);
          end for;
        for  i in 16:18 loop
            C := C + N[i]*delta^d[i]*tao^t[i]*exp(-delta^l[i])*exp(-tao^m[i])*(t[i] -
              m[i]*tao^m[i]);
          end for;
        //P:=rou/(MM)*R*T*(1+A);
        B:=(a[2]+(-0.1)*tao^(-0.1)+a[4]*b[4]/(exp(b[4]*tao)-1)+a[5]*b[5]/(exp(b[5]*tao)-1)+a[6]*b[6]/(exp(b[6]*tao)-1))*tao+1;
        h:=R*T*(B+C+A+1)/MM+6.3121e4;
        end SimpleEnthalpy_Trou;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end R125_USD;
    end MediaUSD;

    package R421mix

       extends MultiPhaseMixture.Templates.ExternalTwoPhaseMixture(
      nP=2,
      setupInfo(
        compounds="",
        libraryName = "RefProp",
        setupInfo="mixture=R421A.mix|path=C:/Program Files (x86)/REFPROP/|eos=default"),
      substanceNames={"R125","R134a"},
      reference_X={0.7,0.3},
      mediumName="r421a");

      extends MultiPhaseMixture.Icons.Media;

      redeclare model extends ThermoProperties
      end ThermoProperties;

     redeclare model extends Properties
     end Properties;

     redeclare model extends MultiPhaseProperties
     end MultiPhaseProperties;

      annotation (Icon(graphics={
          Polygon(
            points={{62,-72},{54,-76},{6,-76},{50,-72},{50,-70},{34,-36},{62,-72}},
            lineColor={255,255,255},
            smooth=Smooth.None,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));

    end R421mix;
  end Media;

  package Components
    package source_sink_ports
      model InFlow
        replaceable package medium =TJUthermo.Media constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model"
        annotation(choicesAllMatching=true);
          parameter Modelica.SIunits.MassFlowRate m0=1;
          parameter Modelica.SIunits.SpecificEnthalpy h0=1.5e5;
          parameter Modelica.SIunits.Pressure p0=1e5;
          parameter Boolean UseT=false;
          parameter Modelica.SIunits.Temperature T;
        Interface.Flange_a node annotation (Placement(transformation(extent={{-86,
                  -14},{-66,6}}), iconTransformation(extent={{-86,-14},{-66,6}})));
      equation
        if UseT then
        node.m_flow=m0;
        node.h=medium.specificEnthalpy_pT(node.p,T);
        else
        node.h=h0;
        node.m_flow=m0;
        //T=medium.temperature_ph(node.p,node.h);
        end if;
       //node.p=p0;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-74,52},{82,-58}},
                lineColor={28,108,200},
                fillColor={0,127,0},
                fillPattern=FillPattern.Solid), Text(
                extent={{-48,-60},{46,-96}},
                lineColor={28,108,200},
                fillColor={0,127,0},
                fillPattern=FillPattern.Solid,
                textString="Inflow")}),                                Diagram(
              coordinateSystem(preserveAspectRatio=false)),
          uses(Modelica(version="3.2.1")));
      end InFlow;

    model OutFlow
    replaceable package medium =TJUthermo.Media constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model"
      annotation(choicesAllMatching=true);
    TJUthermo.Interface.Flange_b node
        annotation (Placement(transformation(extent={{72,-12},{92,8}})));
        parameter Modelica.SIunits.MassFlowRate m0=1;
        parameter Modelica.SIunits.SpecificEnthalpy h0=2.5e5;
        parameter Modelica.SIunits.Pressure p0=1e5;
         Modelica.SIunits.Temperature T;
    equation
    node.p=p0;
     //node.m_flow=m0;
     node.h=h0;
     T=medium.temperature_ph(node.p,inStream(node.h))
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
              extent={{-62,40},{84,-38}},
              lineColor={28,108,200},
              fillColor={0,128,255},
              fillPattern=FillPattern.Solid), Text(
              extent={{-38,-46},{50,-84}},
              lineColor={28,108,200},
              fillColor={0,128,255},
              fillPattern=FillPattern.Solid,
              textString="outflow")}),                               Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        uses(Modelica(version="3.2.1")));
    end OutFlow;

    model Wall
      parameter Modelica.SIunits.Temperature T=350;

      TJUthermo.Interface.ThermoT thermoT
        annotation (Placement(transformation(extent={{-26,-26},{24,24}})));
    equation
      thermoT.T=T;
      annotation (Icon(graphics={Rectangle(
              extent={{-62,30},{68,6}},
              lineColor={28,108,200},
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid), Text(
              extent={{-46,-6},{42,-28}},
              lineColor={28,108,200},
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid,
              textString="Wall")}));
    end Wall;

      model Walls

        Interface.mutiThermoP walls(final N=N)
          annotation (Placement(transformation(extent={{-44,-26},{38,30}})));
          parameter Modelica.SIunits.Temperature T1;
          parameter Modelica.SIunits.Temperature Tn;
          parameter Integer N(min=2)=2;
      equation
        walls.T=linspace(T1,Tn,N);
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-68,20},{62,-18}},
                lineColor={28,108,200},
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid)}), Diagram(coordinateSystem(
                preserveAspectRatio=false)));
      end Walls;

      model DoubleWalls "双面墙"

        Interface.mutiThermoP Wall_side1(final N=N)
          annotation (Placement(transformation(extent={{-62,12},{62,32}})));
        Interface.mutiThermoP Wall_side2(final N=N)
          annotation (Placement(transformation(extent={{-62,-16},{62,6}})));
          parameter Integer N(min=2)=2;
          parameter Modelica.SIunits.Temperature T1;
          parameter Modelica.SIunits.Temperature Tn;
          parameter Modelica.SIunits.Temperature[N] Tstart=linspace(T1,Tn,N);
          Modelica.SIunits.Temperature[N] T;
          parameter Modelica.SIunits.Mass M;
          parameter Modelica.SIunits.SpecificHeatCapacity Cp;
      initial equation
            Wall_side1.T=Tstart;
            Wall_side2.T=Tstart;
      equation
        Wall_side1.T=Wall_side2.T;
        Wall_side1.T=T;
        Wall_side1.Q+Wall_side2.Q=M*Cp*der(T);
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
                extent={{-36,22},{-32,-4}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-22,22},{-18,-4}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-8,22},{-4,-4}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{8,22},{12,-4}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{22,22},{26,-4}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{36,22},{40,-4}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-100,54},{100,-34}},
                lineColor={28,108,200},
                fillColor={170,255,255},
                fillPattern=FillPattern.Solid)}),                      Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
                extent={{-36,22},{-32,-4}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-22,22},{-18,-4}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-8,22},{-4,-4}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{8,22},{12,-4}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{22,22},{26,-4}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{36,22},{40,-4}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-100,50},{100,-38}},
                lineColor={28,108,200},
                fillColor={170,255,255},
                fillPattern=FillPattern.Solid)}));
      end DoubleWalls;

      model InFlow_liquid

        Interface.SimpleFlange_a A
          annotation (Placement(transformation(extent={{84,-12},{104,8}})));
          parameter Modelica.SIunits.Temperature T0;
      equation
        A.T=T0;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-100,48},{100,-56}},
                lineColor={28,108,200},
                fillColor={0,127,127},
                fillPattern=FillPattern.Forward), Text(
                extent={{-70,30},{46,-32}},
                lineColor={213,170,255},
                fillColor={0,127,127},
                fillPattern=FillPattern.Forward,
                textString="Liquid")}), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end InFlow_liquid;

      model InFlow_liquid_useT

        Interface.SimpleFlange_a A
          annotation (Placement(transformation(extent={{84,-12},{104,8}})));
        Modelica.Blocks.Interfaces.RealInput u annotation (Placement(transformation(
                extent={{-30,28},{10,68}}), iconTransformation(
              extent={{-16,-16},{16,16}},
              rotation=-90,
              origin={-12,40})));
      equation
        A.T=u;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-100,48},{100,-56}},
                lineColor={28,108,200},
                fillColor={0,127,127},
                fillPattern=FillPattern.Forward), Text(
                extent={{-66,30},{50,-32}},
                lineColor={213,170,255},
                fillColor={0,127,127},
                fillPattern=FillPattern.Forward,
                textString="InFlow_liquid")}),
                                        Diagram(coordinateSystem(preserveAspectRatio=false)));
      end InFlow_liquid_useT;

      model OutFlow_liquid

        Interface.SimpleFlange_b A
          annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
          parameter Modelica.SIunits.MassFlowRate m0;
         // parameter Modelica.SIunits.Temperature T0;
      equation
        A.m_flow=m0
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-98,54},{96,-56}},
                lineColor={28,108,200},
                fillColor={85,170,255},
                fillPattern=FillPattern.Forward), Text(
                extent={{-66,22},{62,-28}},
                lineColor={215,215,215},
                fillColor={85,170,255},
                fillPattern=FillPattern.Forward,
                textString="liquid")}), Diagram(coordinateSystem(preserveAspectRatio=false),
              graphics={Rectangle(
                extent={{-100,50},{94,-60}},
                lineColor={28,108,200},
                fillColor={0,127,127},
                fillPattern=FillPattern.Forward)}));
       // A.T=T0;
      end OutFlow_liquid;

      model Double_wall

        Interface.ThermoT A
          annotation (Placement(transformation(extent={{-10,14},{10,34}})));
        Interface.ThermoT B
          annotation (Placement(transformation(extent={{-10,-44},{10,-24}})));
          parameter Modelica.SIunits.SpecificHeatCapacity Cp;
          parameter Modelica.SIunits.Temperature T0;
          Modelica.SIunits.Temperature T(start=T0);
          parameter Modelica.SIunits.Mass M;
      equation
        A.T=B.T;
        T=A.T;
        A.Q+B.Q=Cp*M*der(T);
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-100,52},{100,-58}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Forward)}), Diagram(coordinateSystem(
                preserveAspectRatio=false)));
      end Double_wall;

      model CounterDoubleWalls

        Interface.mutiThermoP Wall_side1(final N=N)
          annotation (Placement(transformation(extent={{-62,16},{62,36}}),
              iconTransformation(extent={{-62,16},{62,36}})));
        Interface.mutiThermoP Wall_side2(final N=N)
          annotation (Placement(transformation(extent={{-60,-50},{64,-28}}),
              iconTransformation(extent={{-60,-50},{64,-28}})));
          parameter Integer N(min=2)=2;
          parameter Modelica.SIunits.Temperature T1;
          parameter Modelica.SIunits.Temperature Tn;
          parameter Modelica.SIunits.Temperature[N] Tstart=linspace(T1,Tn,N);
          Modelica.SIunits.Temperature[N] T;
          parameter Modelica.SIunits.Mass M;
          parameter Modelica.SIunits.SpecificHeatCapacity Cp;
          TJUthermo.DataRecord.Arrays summary(n=N,T=T);
      initial equation
            Wall_side1.T=Tstart;
            Wall_side2.T=Tstart[N:-1:1];
      equation
        Wall_side1.T=Wall_side2.T[N:-1:1];
        Wall_side1.T=T;
        Wall_side1.Q+Wall_side2.Q[N:-1:1]=(M/N)*Cp*der(T);
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-98,48},{100,-60}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.HorizontalCylinder), Text(
                extent={{-82,16},{72,-22}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={175,175,175},
                textString="Counter Flow Wall")}),                     Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end CounterDoubleWalls;

      model InFlow_ph
       replaceable package medium =TJUthermo.Media constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model"
        annotation(choicesAllMatching=true);
        Interface.Flange_a node
          annotation (Placement(transformation(extent={{66,-10},{86,10}})));
        parameter Modelica.SIunits.MassFlowRate m0=1;
          parameter Modelica.SIunits.SpecificEnthalpy h0=1.5e5;
          parameter Modelica.SIunits.Pressure p0=1e5;
          parameter Boolean UseT=false;
          parameter Modelica.SIunits.Temperature T;
      equation
        if UseT then
        node.p=p0;
        node.h=medium.specificEnthalpy_pT(node.p,T);
        else
        node.h=h0;
        node.p=p0;
        //T=medium.temperature_ph(node.p,node.h);
        end if;
       //node.p=p0;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-70,44},{78,-52}},
                lineColor={28,108,200},
                fillColor={85,255,255},
                fillPattern=FillPattern.Solid), Text(
                extent={{-56,20},{48,-30}},
                lineColor={28,108,200},
                fillColor={85,255,255},
                fillPattern=FillPattern.Solid,
                textString="inflow_ph")}), Diagram(coordinateSystem(
                preserveAspectRatio=false)));
      end InFlow_ph;

      model OutFlow_mh
      replaceable package medium =TJUthermo.Media constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model"
        annotation(choicesAllMatching=true);
      TJUthermo.Interface.Flange_b node
          annotation (Placement(transformation(extent={{72,-12},{92,8}})));
          parameter Modelica.SIunits.MassFlowRate m0=1;
          parameter Modelica.SIunits.SpecificEnthalpy h0=2.5e5;
          parameter Modelica.SIunits.Pressure p0=1e5;
           Modelica.SIunits.Temperature T;
      equation
      node.p=p0;
       //node.m_flow=m0;
       node.h=h0;
       T=medium.temperature_ph(node.p,inStream(node.h))
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-62,40},{84,-38}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid), Text(
                extent={{-42,22},{46,-16}},
                lineColor={28,108,200},
                fillColor={0,128,255},
                fillPattern=FillPattern.Solid,
                textString="outflow_mh")}),                            Diagram(
              coordinateSystem(preserveAspectRatio=false)),
          uses(Modelica(version="3.2.1")));
      end OutFlow_mh;

      model InFlow_Tmflow

        Interface.SimpleFlange_b outFlow
          annotation (Placement(transformation(extent={{90,-10},{110,10}})));
          parameter Modelica.SIunits.Temperature T0;
          parameter Modelica.SIunits.MassFlowRate m0;
          parameter Boolean use_m0=false;
      equation
          outFlow.T=T0;
          if use_m0 then
          outFlow.m_flow=m0;
          end if;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-100,28},{100,-34}},
                lineColor={28,108,200},
                fillColor={255,85,85},
                fillPattern=FillPattern.Solid), Text(
                extent={{-66,14},{64,-20}},
                lineColor={28,108,200},
                fillColor={255,85,85},
                fillPattern=FillPattern.Solid,
                textString="InFlow_Tmflow")}), Diagram(coordinateSystem(
                preserveAspectRatio=false)));
      end InFlow_Tmflow;

      model OutFlow_Tmflow

        Interface.SimpleFlange_a inFlow
          annotation (Placement(transformation(extent={{-110,-12},{-90,8}})));
          parameter Integer N;
          parameter Modelica.SIunits.MassFlowRate m0;
          parameter Boolean use_m0=true;
      equation
         if use_m0 then
         inFlow.m_flow=m0/N;
         end if;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-100,28},{100,-32}},
                lineColor={28,108,200},
                fillColor={170,213,255},
                fillPattern=FillPattern.Solid), Text(
                extent={{-74,16},{58,-18}},
                lineColor={28,108,200},
                fillColor={170,213,255},
                fillPattern=FillPattern.Solid,
                textString="outFlow_Tmflow")}), Diagram(coordinateSystem(
                preserveAspectRatio=false)));
      end OutFlow_Tmflow;

      model OutFlow_liquid_useM

        Interface.SimpleFlange_b A
          annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
          //parameter Modelica.SIunits.MassFlowRate m0;
         // parameter Modelica.SIunits.Temperature T0;
        Modelica.Blocks.Interfaces.RealInput u annotation (Placement(transformation(
                extent={{-22,32},{18,72}}), iconTransformation(
              extent={{20,-20},{-20,20}},
              rotation=90,
              origin={-2,42})));
      equation
        A.m_flow=u;
       // A.T=T0;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-98,54},{96,-56}},
                lineColor={28,108,200},
                fillColor={85,170,255},
                fillPattern=FillPattern.Forward), Text(
                extent={{-66,22},{62,-28}},
                lineColor={215,215,215},
                fillColor={85,170,255},
                fillPattern=FillPattern.Forward,
                textString="outFlow_liquid")}),
                                        Diagram(coordinateSystem(preserveAspectRatio=false),
              graphics={Rectangle(
                extent={{-100,50},{94,-60}},
                lineColor={28,108,200},
                fillColor={0,127,127},
                fillPattern=FillPattern.Forward)}));
      end OutFlow_liquid_useM;

      model MultiSource_liquid
        parameter Integer N=4;
        Interface.SimpleFlangeMulti outFlow(N=N)
          annotation (Placement(transformation(extent={{90,-10},{110,10}})));

          parameter Modelica.SIunits.Temperature T0;
          parameter Modelica.SIunits.MassFlowRate m0;
          parameter Boolean use_m0=false;
      equation
        for i in 1:N loop
          outFlow.T[i]=T0;
        end for;
          if use_m0 then
            for i in 1:N loop
          outFlow.m_flow[i]=m0/N;
            end for;
          end if;
          annotation (Placement(transformation(extent={{58,-12},{78,8}})),
                    Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-82,42},{98,-60}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid), Text(
                extent={{-48,6},{50,-40}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid,
                textString="Source")}),                                Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
                extent={{-78,38},{68,-46}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid), Text(
                extent={{-60,8},{36,-26}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid,
                textString="InFlow")}));
      end MultiSource_liquid;

      model MultiSource_liquid_sink
        parameter Integer N=4;
          parameter Modelica.SIunits.MassFlowRate m0;
          parameter Boolean use_m0=true;
        Interface.SimpleFlangeMulti inFlow(N=N)
          annotation (Placement(transformation(extent={{-106,-10},{-86,10}}),
              iconTransformation(extent={{-106,-10},{-86,10}})));
      equation
         if use_m0 then
           for i in 1:N loop
         inFlow.m_flow[i]=m0;
           end for;
         end if;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-92,30},{92,-46}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid), Text(
                extent={{-74,10},{46,-26}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid,
                textString="Sink")}),                                  Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end MultiSource_liquid_sink;

      model R134AInFlow
        replaceable package medium=ExternalMedia.Media.CoolPropMedium (mediumName="R134A");
        TJUthermo.Interface.Flange_a node
          annotation (Placement(transformation(extent={{72,-12},{92,8}})));
          parameter Modelica.SIunits.MassFlowRate m0=1;
          parameter Modelica.SIunits.SpecificEnthalpy h0=1.5e5;
          parameter Modelica.SIunits.Pressure p0=1e5;
          parameter Boolean UseT=false;
          parameter Modelica.SIunits.Temperature T;
      equation
        if UseT then
        node.m_flow=m0;
        node.h=medium.specificEnthalpy_pT(node.p,T);
        else
        node.h=h0;
        node.m_flow=m0;
        //T=medium.temperature_ph(node.p,node.h);
        end if;
       //node.p=p0;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-74,52},{82,-58}},
                lineColor={28,108,200},
                fillColor={0,127,0},
                fillPattern=FillPattern.Solid), Text(
                extent={{-48,-60},{46,-96}},
                lineColor={28,108,200},
                fillColor={0,127,0},
                fillPattern=FillPattern.Solid,
                textString="Inflow")}),                                Diagram(
              coordinateSystem(preserveAspectRatio=false)),
          uses(Modelica(version="3.2.1")));
      end R134AInFlow;

      model R134AOutFlow
      replaceable package medium=ExternalMedia.Media.CoolPropMedium (mediumName="R134A");
      TJUthermo.Interface.Flange_b node
          annotation (Placement(transformation(extent={{72,-12},{92,8}})));
          parameter Modelica.SIunits.MassFlowRate m0=1;
          parameter Modelica.SIunits.SpecificEnthalpy h0=2.5e5;
          parameter Modelica.SIunits.Pressure p0=1e5;
           Modelica.SIunits.Temperature T;
      equation
      node.p=p0;
       //node.m_flow=m0;
       node.h=h0;
       T=medium.temperature_ph(node.p,inStream(node.h))
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-62,40},{84,-38}},
                lineColor={28,108,200},
                fillColor={0,128,255},
                fillPattern=FillPattern.Solid), Text(
                extent={{-38,-46},{50,-84}},
                lineColor={28,108,200},
                fillColor={0,128,255},
                fillPattern=FillPattern.Solid,
                textString="outflow")}),                               Diagram(
              coordinateSystem(preserveAspectRatio=false)),
          uses(Modelica(version="3.2.1")));

      end R134AOutFlow;
    end source_sink_ports;

    package units
      model Cell
      replaceable package medium =TJUthermo.Media constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model"
        annotation(choicesAllMatching=true);
        TJUthermo.Interface.Flange_a A
          annotation (Placement(transformation(extent={{-108,-10},{-88,10}})));
        TJUthermo.Interface.Flange_b B
          annotation (Placement(transformation(extent={{88,-10},{108,10}})));
        TJUthermo.Interface.ThermoT wall
          annotation (Placement(transformation(extent={{-14,8},{12,34}})));
          //常量区//
          constant Real pi=3.1415;
             parameter Modelica.SIunits.CoefficientOfHeatTransfer Ul=1500
          "液态换热系数";
             parameter Modelica.SIunits.CoefficientOfHeatTransfer Utp=5000
          "两相换热系数";
             parameter Modelica.SIunits.CoefficientOfHeatTransfer Uv=1000
          "气态换热系数";
             parameter Modelica.SIunits.Area F=0.0063;
             parameter Modelica.SIunits.Volume V=0.7263;
             parameter Real max_drdp=1;
             parameter Real max_drdh=1;
             parameter Real witdth_x=0.05 "减缓突变设定值";
             parameter Modelica.SIunits.SpecificEnthalpy hin_start=medium.specificEnthalpy_pT(pstart,Tin_start) annotation(Dialog(tab="Initialization"));
             parameter Modelica.SIunits.SpecificEnthalpy hout_start=medium.specificEnthalpy_pT(pstart,Tout_start) annotation(Dialog(tab="Initialization"));
            // parameter Modelica.SIunits.Pressure pstart=1.2e6;
            parameter Modelica.SIunits.MassFlowRate m_flow=0.5;
            parameter Modelica.SIunits.Temperature Tin_start=313 annotation(Dialog(tab="Initialization"));
            parameter Modelica.SIunits.Temperature Tout_start=320 annotation(Dialog(tab="Initialization"));
            parameter Modelica.SIunits.Pressure pstart=1.2e6 annotation(Dialog(tab="Initialization"));
            parameter Boolean Mcons=false;
                  //parameter Modelica.SIunits.Mass M=94;
          //变量区//
          Modelica.SIunits.CoefficientOfHeatTransfer U;
          Modelica.SIunits.SpecificEnthalpy h;//(start=0.5*(hin_start+hout_start));
          Modelica.SIunits.SpecificEnthalpy hin;//(start=hin_start);
          Modelica.SIunits.SpecificEnthalpy hl "饱和液焓值";
          Modelica.SIunits.SpecificEnthalpy hv "饱和气焓值";
          Real x "名义干度";
          Modelica.SIunits.SpecificEnthalpy hout;//(start=hout_start);
          Modelica.SIunits.MassFlowRate m_in;//(start=m_flow);
          Modelica.SIunits.MassFlowRate m_out;//(start=-m_flow);
          Modelica.SIunits.Pressure p(start=pstart);
          Modelica.SIunits.HeatFlowRate Q;//【w】
          Modelica.SIunits.HeatFlowRate Q_final;//【w】
          Modelica.SIunits.Temperature Twall "K";
          Modelica.SIunits.Temperature T(start=0.5*(Tin_start+Tout_start));
          medium.SaturationProperties Sat1;
          medium.ThermodynamicState state;
          Modelica.SIunits.Density rou;//(start=medium.density_ph(pstart,0.5*(hin_start+hout_start))) "main";
          Modelica.SIunits.DerDensityByEnthalpy drdh;
          Modelica.SIunits.DerDensityByPressure drdp;
          Modelica.SIunits.Mass M;

      initial equation
        h=0.5*(hin_start+hout_start);
        rou=medium.density_ph(pstart,0.5*(hin_start+hout_start)) "main";
        //rou=state.d;
       // h=hstart;
       // T=Tstart;
      equation
        U=TJUthermo.HeatTransferModel.constHeatTransfer(Utp=Utp,Ul=Ul,Uv=Uv,x=x,witdth_x=witdth_x);
        Sat1=medium.setSat_p(p);
        hl=Sat1.hl;
        hv=Sat1.hv;
        x=(h-hl)/(hv-hl);
        M=rou*V;
        drdh=max(max_drdh/(-4000),medium.density_derh_p(state));
        drdp=min(max_drdp/1e5,medium.density_derp_h(state));
        state = medium.setState_ph(p,h);
        T=state.T;
        Twall=wall.T;
        wall.Q=Q;
        Q_final=m_out*hout+m_in*hin;
        Q=U*F*(Twall-T);
        m_in=A.m_flow;
        m_out=B.m_flow;
        if Mcons then
          m_in+m_out=0;
        else
        m_in+m_out=V*der(rou);//======质量平衡=======//
        end if;
        der(rou)=drdp*der(p)+drdh*der(h);
          A.h=inStream(A.h);
          B.h=hout;
          m_in*hin+m_out*hout+Q=M*der(h)-V*der(p);//2*h=hin+hout;能量平衡
          hin=inStream(A.h);
          hout=2*h-hin;
          A.p=p;
          B.p=p;
        //  p=1e5;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-92,34},{94,-40}},
                lineColor={0,0,0},
                fillColor={170,170,255},
                fillPattern=FillPattern.Solid), Text(
                extent={{-58,-52},{50,-76}},
                lineColor={0,0,0},
                fillColor={170,170,255},
                fillPattern=FillPattern.Solid,
                textString="换热管道")}), Diagram(coordinateSystem(
                preserveAspectRatio=false)));
      end Cell;

     model SimpleFlow
      replaceable package Medium =TJUthermo.Media constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model"
       annotation(choicesAllMatching=true);
       //接口//
      TJUthermo.Interface.Flange_a InFlow
        annotation (Placement(transformation(extent={{-106,-6},{-86,14}})));
      TJUthermo.Interface.Flange_b OutFlow
        annotation (Placement(transformation(extent={{88,-4},{108,16}})));

      TJUthermo.Interface.mutiThermoP Wall(N=N)
        annotation (Placement(transformation(extent={{-34,26},{38,70}})));
      //数据汇总//
      TJUthermo.DataRecord.SummaryClass Summary(T_profiles(n=N,T=Cell[:].T),N=N,h_cell=Cell[:].h,h_node=h_node_,M_node=M_node_,x=Cell.x,p=Cell[1].p);
      //常量//
      Modelica.SIunits.HeatFlowRate Q;
      parameter Integer N(min=1)=2 "分割单元体数目";
          parameter Modelica.SIunits.CoefficientOfHeatTransfer Ul=1000
          "液态换热系数";
           parameter Modelica.SIunits.CoefficientOfHeatTransfer Utp=5000
          "两相换热系数";
           parameter Modelica.SIunits.CoefficientOfHeatTransfer Uv=2000
          "气态换热系数";
               parameter Boolean Mcons=false;
           parameter Modelica.SIunits.Area F=0.0063;
           parameter Modelica.SIunits.Volume V=0.7263;
           parameter Real max_drdp=1;
           parameter Real max_drdh=1;
           parameter Real witdth_x=0.05 "减缓突变设定值";
           //初值设定区//
           parameter Modelica.SIunits.Temperature Tin_start=303 annotation(Dialog(tab="初始化参数"));
           parameter Modelica.SIunits.Temperature Tout_start annotation(Dialog(tab="初始化参数"));
           parameter Modelica.SIunits.MassFlowRate m_flow_start=0.5 annotation(Dialog(tab="初始化参数"));
           parameter Modelica.SIunits.Temperature Tstart[N+1]=Medium.temperature_ph(pstart,hstart);
           parameter Modelica.SIunits.SpecificEnthalpy hstart[N+1]=linspace(Medium.specificEnthalpy_pT(pstart,Tin_start),Medium.specificEnthalpy_pT(pstart,Tout_start),N+1);
           parameter Modelica.SIunits.Pressure pstart=1.2e6 annotation(Dialog(tab="初始化参数"));
           //parameter Modelica.SIunits.Temperature Tstart=313;
          //变量//
       TJUthermo.Interface.ThermoConvet ThermoC(final N=N);
        replaceable TJUthermo.Components.units.Cell Cell[N](
        redeclare each final package medium=Medium,
        each Mcons=Mcons,
        each Ul=Ul,
        each Utp=Utp,
        each Uv=Uv,
        each F=F/N,
        each V=V/N,
        each m_flow=m_flow_start,
        final Tin_start=Tstart[1:N],
        final Tout_start=Tstart[2:N+1],
        final hin_start=hstart[1:N],
        final hout_start=hstart[2:N+1],
        each pstart=pstart,
        each max_drdp=max_drdp,
        each max_drdh=max_drdh,
        each witdth_x=witdth_x);
        //Tstart=Tstart,
      protected
        Modelica.SIunits.SpecificEnthalpy[N+1] h_node_;
        Modelica.SIunits.MassFlowRate[N+1] M_node_;
       // Real[N] x;
     equation
       //单元体首尾连接//
      connect(InFlow,Cell[1].A);
      for i in 1:N-1 loop
       connect(Cell[i].B,Cell[i+1].A);
      end for;
      connect(Cell[N].B,OutFlow);
      //单元体壁面连接//
      connect(Cell.wall,ThermoC.single);
      connect(Wall,ThermoC.multi);
      //数据汇总//
      h_node_[1:N]=Cell[:].hin;
      h_node_[N+1]=Cell[N].hout;
      M_node_[1:N]=Cell[:].m_in;
      M_node_[N+1]=-Cell[N].m_out;//change sign
     algorithm
       Q:=0;
       for i in 1:N loop
        Q:=Q + Cell[i].Q_final;
       end for;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
              extent={{-90,42},{92,-28}},
              lineColor={28,108,200},
              fillColor={170,255,255},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{-100,42},{-66,-28}},
              lineColor={28,108,200},
              fillColor={0,128,255},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{66,42},{100,-28}},
              lineColor={28,108,200},
              fillColor={0,128,255},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{-72,16},{40,-2}},
              lineColor={28,108,200},
              fillColor={0,127,0},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{20,26},{46,6},{22,-18},{66,6},{20,26}},
              lineColor={28,108,200},
              fillColor={0,127,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-56,40},{14,20}},
              lineColor={28,108,200},
              fillColor={0,127,0},
              fillPattern=FillPattern.Solid,
              textString="Flow")}),                                  Diagram(
            coordinateSystem(preserveAspectRatio=false)));
     end SimpleFlow;

      model Cell_liquid

        TJUthermo.Interface.SimpleFlange_a A
          annotation (Placement(transformation(extent={{-108,-10},{-88,10}})));
        TJUthermo.Interface.SimpleFlange_b B
          annotation (Placement(transformation(extent={{88,-10},{108,10}})));
        TJUthermo.Interface.ThermoT wall
          annotation (Placement(transformation(extent={{-14,8},{12,34}})));
          //常量区//
          constant Real pi=3.1415;
          parameter Modelica.SIunits.Area F=0.0063;
          parameter Modelica.SIunits.Volume V=0.7263;
          Modelica.SIunits.MassFlowRate m_flow;//修改处
          parameter Modelica.SIunits.Temperature Tstart=313 annotation(Dialog(tab="Initialization"));
          parameter Modelica.SIunits.CoefficientOfHeatTransfer  U;
          parameter Modelica.SIunits.SpecificHeatCapacity Cp;
          parameter Modelica.SIunits.Mass M;
          Modelica.SIunits.HeatFlowRate Q;//【w】
          Modelica.SIunits.Temperature Twall "K";
          Modelica.SIunits.Temperature T;//(start=Tstart);
          Modelica.SIunits.Temperature Tin;//(start=Tstart);
          Modelica.SIunits.Temperature Tout;//(start=Tstart);
      initial equation
        T=Tstart;

      equation
        Twall=wall.T;
        wall.Q=Q;
        Q=U*F*(Twall-T);
        m_flow=A.m_flow;
        A.m_flow+B.m_flow=0;
        Tin=A.T;
        B.T=Tout;
        m_flow*Cp*(Tin-Tout)+Q=M*Cp*der(T);//2*h=hin+hout;
        Tout=2*T-Tin;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-96,38},{96,-50}},
                lineColor={28,108,200},
                fillColor={85,170,255},
                fillPattern=FillPattern.Forward), Text(
                extent={{-58,6},{50,-38}},
                lineColor={28,108,200},
                fillColor={85,170,255},
                fillPattern=FillPattern.Forward,
                textString="Liquid")}),                                Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end Cell_liquid;

      model Cell_liquid_type2

        TJUthermo.Interface.SimpleFlange_a A
          annotation (Placement(transformation(extent={{-108,-10},{-88,10}})));
        TJUthermo.Interface.SimpleFlange_b B
          annotation (Placement(transformation(extent={{88,-10},{108,10}})));
        TJUthermo.Interface.ThermoT wall
          annotation (Placement(transformation(extent={{-14,8},{12,34}})));
          //常量区//
          constant Real pi=3.1415;
          parameter Modelica.SIunits.Area F=0.0063;
          parameter Modelica.SIunits.Volume V=0.7263;
          Modelica.SIunits.MassFlowRate m_flow;//修改处
          parameter Modelica.SIunits.Temperature Tstart=313 annotation(Dialog(tab="Initialization"));
          parameter Modelica.SIunits.CoefficientOfHeatTransfer  U;
          parameter Modelica.SIunits.SpecificHeatCapacity Cp;
          parameter Modelica.SIunits.Mass M;
          Modelica.SIunits.HeatFlowRate Q;//【w】
          Modelica.SIunits.Temperature Twall "K";
          Modelica.SIunits.Temperature T(start=Tstart);
          Modelica.SIunits.Temperature Tin;//(start=Tstart);
          Modelica.SIunits.Temperature Tout;//(start=Tstart);
      //initial equation
       // T=Tstart;

      equation
        Twall=wall.T;
        wall.Q=Q;
        Q=U*F*(Twall-T);
        m_flow=A.m_flow;
        A.m_flow+B.m_flow=0;
        Tin=A.T;
        B.T=Tout;
        m_flow*Cp*(Tin-Tout)+Q=M*Cp*der(T);//2*h=hin+hout;
        Tout=2*T-Tin;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-96,38},{96,-50}},
                lineColor={28,108,200},
                fillColor={85,170,255},
                fillPattern=FillPattern.Forward), Text(
                extent={{-58,6},{50,-38}},
                lineColor={28,108,200},
                fillColor={85,170,255},
                fillPattern=FillPattern.Forward,
                textString="Liquid")}),                                Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end Cell_liquid_type2;

      model SimpleFlowliquid

        Interface.SimpleFlange_a node_a
          annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
        Interface.SimpleFlange_b node_b
          annotation (Placement(transformation(extent={{88,-8},{108,12}})));
        Interface.mutiThermoP node_wall(N=N)
          annotation (Placement(transformation(extent={{-46,10},{40,52}})));
          //数据汇总//
        TJUthermo.DataRecord.Arrays Summary(n=N,T=Cell[:].T);
          //常量//
          parameter Integer N(min=2)=2;
          parameter Modelica.SIunits.Area F=0.0063;
          parameter Modelica.SIunits.Volume V=0.7263;
          Modelica.SIunits.MassFlowRate m_flow=0.5;
          parameter Modelica.SIunits.Temperature T1=313 annotation(Dialog(tab="Initialization"));
          parameter Modelica.SIunits.Temperature Tn=313 annotation(Dialog(tab="Initialization"));
          parameter Modelica.SIunits.Temperature Tstart[N]=linspace(T1,Tn,N) annotation(Dialog(tab="Initialization"));
          parameter Modelica.SIunits.CoefficientOfHeatTransfer  U;
          parameter Modelica.SIunits.SpecificHeatCapacity Cp;
          parameter Modelica.SIunits.Mass M;

          //变量//
          Modelica.SIunits.HeatFlowRate Q;
          replaceable TJUthermo.Components.units.Cell_liquid Cell[N](
          each F=F/N,
          each V=V/N,
          each U=U,
          each Cp=Cp,
          each M=M/N,
          Tstart=Tstart);

          TJUthermo.Interface.ThermoConvet ThermoC(final N=N);
      equation
        connect(node_a,Cell[1].A);
        for i in 1:N-1 loop
          connect(Cell[i].B,Cell[i+1].A);
        end for;
        connect(Cell[N].B,node_b);
        //连接壁面温度//
        connect(Cell.wall,ThermoC.single);
        connect(ThermoC.multi,node_wall);
      algorithm
        Q:=0;
        for i in 1:N loop
        Q:=Q + Cell[i].Q;
        end for;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-100,46},{100,-42}},
                lineColor={28,108,200},
                fillColor={170,170,255},
                fillPattern=FillPattern.Forward), Text(
                extent={{-68,10},{62,-30}},
                lineColor={28,108,200},
                fillColor={170,170,255},
                fillPattern=FillPattern.Solid,
                textString="LiquidFlow"),
              Rectangle(
                extent={{-86,18},{62,10}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={175,175,175}),
              Polygon(
                points={{58,28},{84,16},{60,2},{60,12},{58,28}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={175,175,175})}), Diagram(coordinateSystem(
                preserveAspectRatio=false)));
      end SimpleFlowliquid;

      model CellChoiceMedia
        replaceable package medium =TJUthermo.Media2Type
        constrainedby Modelica.Media.Interfaces.PartialMedium "Medium model"
        annotation(choicesAllMatching=true);
        //constrainedby
        //Modelica.Media.Interfaces.PartialMedium "Medium model"
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-50,40},{52,-40}},
                lineColor={28,108,200},
                fillPattern=FillPattern.VerticalCylinder,
                fillColor={170,170,255}), Ellipse(
                extent={{40,40},{-40,-40}},
                lineColor={28,108,200},
                fillPattern=FillPattern.VerticalCylinder,
                fillColor={175,175,175})}),                            Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end CellChoiceMedia;

      model Cell_collector

        TJUthermo.Interface.SimpleFlange_a inFlow
          annotation (Placement(transformation(extent={{-108,-6},{-88,14}})));
        TJUthermo.Interface.SimpleFlange_b outFlow
          annotation (Placement(transformation(extent={{88,-6},{108,14}})));
          parameter Modelica.SIunits.Area A;
          parameter Real K;
          parameter Modelica.SIunits.Mass M;
          parameter Modelica.SIunits.SpecificHeatCapacity Cp;
          Modelica.SIunits.MassFlowRate m_flow;
          parameter Modelica.SIunits.Temperature T0 "初始温度";
          Modelica.SIunits.Temperature dt;
          Modelica.SIunits.Temperature Tin;
          Modelica.SIunits.Temperature Tout;
          Modelica.SIunits.Temperature T(start=T0);
          Modelica.SIunits.Temperature Tam;
          Modelica.SIunits.HeatFlowRate Q;//【w】
          Real eta;
        Modelica.SIunits.HeatFlux I;
                 //输入变量数//
        Modelica.Blocks.Interfaces.RealInput windspeed annotation (Placement(transformation(
              extent={{-17,-17},{17,17}},
              rotation=270,
              origin={-65,33}), iconTransformation(
              extent={{-9,-9},{9,9}},
              rotation=270,
              origin={-57,25})));
        Modelica.Blocks.Interfaces.RealInput I_in annotation (Placement(transformation(
              extent={{-20,-20},{20,20}},
              rotation=-90,
              origin={46,36}), iconTransformation(
              extent={{-9,-9},{9,9}},
              rotation=-90,
              origin={57,25})));
        Modelica.Blocks.Interfaces.RealInput Tam_in annotation (Placement(transformation(
                extent={{-28,22},{12,62}}), iconTransformation(
              extent={{-8,-8},{8,8}},
              rotation=-90,
              origin={-2,28})));
      //initial equation
        //T=T0;
      equation
          if cardinality(windspeed) == 0 then
          windspeed = 3 "风速";
          end if;
          if cardinality(Tam_in) == 0 then
          Tam_in = 25+273.15 "外界温度";
          end if;
          if cardinality(I_in) == 0 then
          I_in = 700 "外界辐射量";
          end if;
       eta=(K*(73.6-0.004206*dt)+7.44*(dt/I)-0.00958*(dt*dt/I))*0.01;
       Tin=inFlow.T;
       Tout=outFlow.T;
       Tam=Tam_in+273.15;
       dt=0.5*(Tin+Tout)-Tam;
       I=I_in;
       Q=K*A*I*eta;
       Q+m_flow*Cp*(Tin-Tout)=M*Cp*der(T);
       m_flow=inFlow.m_flow;
       inFlow.m_flow+outFlow.m_flow=0;
       Tin+Tout=2*T;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-98,28},{98,-22}},
                lineColor={28,108,200},
                fillColor={255,255,170},
                fillPattern=FillPattern.Solid), Text(
                extent={{-80,18},{78,-12}},
                lineColor={28,108,200},
                fillColor={255,255,170},
                fillPattern=FillPattern.Solid,
                textString="collector unit"),
              Text(
                extent={{-74,22},{-60,18}},
                lineColor={28,108,200},
                fillColor={255,255,170},
                fillPattern=FillPattern.Solid,
                textString="风速"),
              Text(
                extent={{42,22},{52,16}},
                lineColor={28,108,200},
                fillColor={255,255,170},
                fillPattern=FillPattern.Solid,
                textString="辐射量"),
              Text(
                extent={{8,24},{26,16}},
                lineColor={28,108,200},
                fillColor={255,255,170},
                fillPattern=FillPattern.Solid,
                textString="外界温度")}),   Diagram(coordinateSystem(
                preserveAspectRatio=false)));
      end Cell_collector;

      model SimpleFlowliquid_type2

        Interface.SimpleFlange_a node_a
          annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
        Interface.SimpleFlange_b node_b
          annotation (Placement(transformation(extent={{88,-8},{108,12}})));
        Interface.mutiThermoP node_wall(N=N)
          annotation (Placement(transformation(extent={{-46,10},{40,52}})));
          //数据汇总//
        TJUthermo.DataRecord.Arrays Summary(n=N,T=Cell[:].T);
          //常量//
          parameter Integer N(min=2)=2;
          parameter Modelica.SIunits.Area F=0.0063;
          parameter Modelica.SIunits.Volume V=0.7263;
          Modelica.SIunits.MassFlowRate m_flow=0.5;
          parameter Modelica.SIunits.Temperature T1=313 annotation(Dialog(tab="Initialization"));
          parameter Modelica.SIunits.Temperature Tn=313 annotation(Dialog(tab="Initialization"));
          parameter Modelica.SIunits.Temperature Tstart[N]=linspace(T1,Tn,N) annotation(Dialog(tab="Initialization"));
          parameter Modelica.SIunits.CoefficientOfHeatTransfer  U;
          parameter Modelica.SIunits.SpecificHeatCapacity Cp;
          parameter Modelica.SIunits.Mass M;

          //变量//
          Modelica.SIunits.HeatFlowRate Q;
          replaceable TJUthermo.Components.units.Cell_liquid_type2 Cell[N](
          each F=F/N,
          each V=V/N,
          each U=U,
          each Cp=Cp,
          each M=M/N,
          Tstart=Tstart);

          TJUthermo.Interface.ThermoConvet ThermoC(final N=N);
      equation
        connect(node_a,Cell[1].A);
        for i in 1:N-1 loop
          connect(Cell[i].B,Cell[i+1].A);
        end for;
        connect(Cell[N].B,node_b);
        //连接壁面温度//
        connect(Cell.wall,ThermoC.single);
        connect(ThermoC.multi,node_wall);
      algorithm
        Q:=0;
        for i in 1:N loop
        Q:=Q + Cell[i].Q;
        end for;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-100,46},{100,-42}},
                lineColor={28,108,200},
                fillColor={170,170,255},
                fillPattern=FillPattern.Forward), Text(
                extent={{-68,10},{62,-30}},
                lineColor={28,108,200},
                fillColor={170,170,255},
                fillPattern=FillPattern.Solid,
                textString="LiquidFlow2"),
              Rectangle(
                extent={{-86,18},{62,10}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={175,175,175}),
              Polygon(
                points={{58,28},{84,16},{60,2},{60,12},{58,28}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={175,175,175})}), Diagram(coordinateSystem(
                preserveAspectRatio=false)));
      end SimpleFlowliquid_type2;

      model Cell_varyHTC
        replaceable package medium=ExternalMedia.Media.CoolPropMedium (mediumName="R134A");
        TJUthermo.Interface.Flange_a A
          annotation (Placement(transformation(extent={{-108,-10},{-88,10}})));
        TJUthermo.Interface.Flange_b B
          annotation (Placement(transformation(extent={{88,-10},{108,10}})));
        TJUthermo.Interface.ThermoT wall
          annotation (Placement(transformation(extent={{-14,8},{12,34}})));
          //常量区//
          constant Real pi=3.1415;
          parameter Integer np(min=1)=1;
             parameter Modelica.SIunits.Length L=1;
             parameter Modelica.SIunits.Length dh=0.03;
             parameter Real max_drdp=1;
             parameter Real max_drdh=1;
             parameter Real witdth_x=0.05 "减缓突变设定值";
             parameter Modelica.SIunits.SpecificEnthalpy hin_start=medium.specificEnthalpy_pT(pstart,Tin_start) annotation(Dialog(tab="Initialization"));
             parameter Modelica.SIunits.SpecificEnthalpy hout_start=medium.specificEnthalpy_pT(pstart,Tout_start) annotation(Dialog(tab="Initialization"));
            // parameter Modelica.SIunits.Pressure pstart=1.2e6;
            parameter Modelica.SIunits.MassFlowRate m_flow=0.5;
            parameter Modelica.SIunits.Temperature Tin_start=313 annotation(Dialog(tab="Initialization"));
            parameter Modelica.SIunits.Temperature Tout_start=320 annotation(Dialog(tab="Initialization"));
            parameter Modelica.SIunits.Pressure pstart=1.2e6 annotation(Dialog(tab="Initialization"));
            parameter Boolean Mcons=false;
            //parameter Modelica.SIunits.Mass M=94;
          //变量区//
          Modelica.SIunits.Area F;
          Modelica.SIunits.Volume V;
          Modelica.SIunits.CoefficientOfHeatTransfer U;
          Modelica.SIunits.CoefficientOfHeatTransfer U1;
          Modelica.SIunits.SpecificEnthalpy h;//(start=0.5*(hin_start+hout_start));
          Modelica.SIunits.SpecificEnthalpy hin;//(start=hin_start);
          Modelica.SIunits.SpecificEnthalpy hl "饱和液焓值";
          Modelica.SIunits.SpecificEnthalpy hv "饱和气焓值";
          Real x "名义干度";
          Modelica.SIunits.SpecificEnthalpy hout;//(start=hout_start);
          Modelica.SIunits.MassFlowRate m_in;//(start=m_flow);
          Modelica.SIunits.MassFlowRate m_out;//(start=-m_flow);
          Modelica.SIunits.Pressure p(start=pstart);
          Modelica.SIunits.Pressure Pc;
          Modelica.SIunits.HeatFlowRate Q;//【w】
          Modelica.SIunits.Temperature Twall "K";
          Modelica.SIunits.Temperature T(start=0.5*(Tin_start+Tout_start));
          medium.SaturationProperties Sat1;
          medium.ThermodynamicState state;
          Modelica.SIunits.Density rou;//(start=medium.density_ph(pstart,0.5*(hin_start+hout_start))) "main";
          Modelica.SIunits.DerDensityByEnthalpy drdh;
          Modelica.SIunits.DerDensityByPressure drdp;
          Modelica.SIunits.Mass M;

      initial equation
        h=0.5*(hin_start+hout_start);
        rou=medium.density_ph(pstart,0.5*(hin_start+hout_start)) "main";
        //rou=state.d;
       // h=hstart;
       // T=Tstart;
      equation
        Pc=medium.getCriticalPressure();
        F=dh*pi*L;
        V=dh*dh*pi/4*L;
        U1=HeatTransferModel.GnielinskiFcn_medium(mor=m_in,dh=dh,T=Sat1.Tsat-0.5,P=p);
        if h<hl or h>hv then
        U=HeatTransferModel.GnielinskiFcn_medium(mor=m_in,dh=dh,T=T,P=p);
       // U=HeatTransferModel.GnielinskiFcn(mor=mor,dh=dh,V=V,Prf=Prf,lamda=lamda);
        else
        U=HeatTransferModel.ShahFcn(x=0.5,P=p,Pc=Pc,Ul=U1);
        end if;
        //U=TJUthermo.HeatTransferModel.constHeatTransfer(Utp=Utp,Ul=Ul,Uv=Uv,x=x,witdth_x=witdth_x);
        Sat1=medium.setSat_p(p);
        hl=Sat1.hl;
        hv=Sat1.hv;
        x=(h-hl)/(hv-hl);
        M=rou*V;
        drdh=max(max_drdh/(-4000),medium.density_derh_p(state));
        drdp=min(max_drdp/1e5,medium.density_derp_h(state));
        state = medium.setState_ph(p,h);
        T=state.T;
        Twall=wall.T;
        wall.Q=Q;

        Q=U*F*(Twall-T);
        m_in=A.m_flow/np;
        m_out=B.m_flow/np;
        if Mcons then
          m_in+m_out=0;
        else
        m_in+m_out=V*der(rou);//======质量平衡=======//
        end if;
        der(rou)=drdp*der(p)+drdh*der(h);
          A.h=inStream(A.h);
          B.h=hout;
          m_in*hin+m_out*hout+Q=M*der(h)-V*der(p);//2*h=hin+hout;能量平衡
          hin=inStream(A.h);
          hout=2*h-hin;
          A.p=p;
          B.p=p;
        //  p=1e5;
        annotation(choicesAllMatching=true,
                    Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-92,34},{94,-40}},
                lineColor={0,0,0},
                fillColor={170,170,255},
                fillPattern=FillPattern.Solid), Text(
                extent={{-58,-52},{50,-76}},
                lineColor={0,0,0},
                fillColor={170,170,255},
                fillPattern=FillPattern.Solid,
                textString="换热管道"),
              Text(
                extent={{-56,12},{52,-26}},
                lineColor={28,108,200},
                textString="Vary")}),         Diagram(coordinateSystem(
                preserveAspectRatio=false)));
      end Cell_varyHTC;

      model veryFlow
      replaceable package Medium=ExternalMedia.Media.CoolPropMedium (mediumName="R134A");
        //接口//
       TJUthermo.Interface.Flange_a InFlow
         annotation (Placement(transformation(extent={{-106,-6},{-86,14}})));
       TJUthermo.Interface.Flange_b OutFlow
         annotation (Placement(transformation(extent={{88,-4},{108,16}})));

       TJUthermo.Interface.mutiThermoP Wall(N=N)
         annotation (Placement(transformation(extent={{-34,26},{38,70}})));
       //数据汇总//
       TJUthermo.DataRecord.SummaryClass Summary(T_profiles(n=N,T=Cell[:].T),N=N,h_cell=Cell[:].h,h_node=h_node_,M_node=M_node_,x=Cell.x,p=Cell[1].p);
       //常量//
       Modelica.SIunits.HeatFlowRate Q;
       parameter Integer N(min=1)=2 "分割单元体数目";
           parameter Integer np(min=1)=1;
            parameter Modelica.SIunits.Length L=1;
            parameter Modelica.SIunits.Length dh=0.03;
            parameter Real max_drdp=1;
            parameter Real max_drdh=1;
            parameter Real witdth_x=0.05 "减缓突变设定值";
            //初值设定区//
            parameter Modelica.SIunits.Temperature Tin_start=303 annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.Temperature Tout_start annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.MassFlowRate m_flow_start=0.5 annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.Temperature Tstart[N+1]=Medium.temperature_ph(pstart,hstart);
            parameter Modelica.SIunits.SpecificEnthalpy hstart[N+1]=linspace(Medium.specificEnthalpy_pT(pstart,Tin_start),Medium.specificEnthalpy_pT(pstart,Tout_start),N+1);
            parameter Modelica.SIunits.Pressure pstart=1.2e6 annotation(Dialog(tab="初始化参数"));
            //parameter Modelica.SIunits.Temperature Tstart=313;
           //变量//
        TJUthermo.Interface.ThermoConvet ThermoC(final N=N);
         replaceable TJUthermo.Components.units.Cell_varyHTC Cell[N](
         redeclare each final package medium=Medium,
         each L=L/N,
         each np=np,
         each dh=dh,
         each m_flow=m_flow_start,
         final Tin_start=Tstart[1:N],
         final Tout_start=Tstart[2:N+1],
         final hin_start=hstart[1:N],
         final hout_start=hstart[2:N+1],
         each pstart=pstart,
         each max_drdp=max_drdp,
         each max_drdh=max_drdh,
         each witdth_x=witdth_x);
         //Tstart=Tstart,
      protected
         Modelica.SIunits.SpecificEnthalpy[N+1] h_node_;
         Modelica.SIunits.MassFlowRate[N+1] M_node_;
        // Real[N] x;
      equation
        //单元体首尾连接//
       connect(InFlow,Cell[1].A);
       for i in 1:N-1 loop
        connect(Cell[i].B,Cell[i+1].A);
       end for;
       connect(Cell[N].B,OutFlow);
       //单元体壁面连接//
       connect(Cell.wall,ThermoC.single);
       connect(Wall,ThermoC.multi);
       //数据汇总//
       h_node_[1:N]=Cell[:].hin;
       h_node_[N+1]=Cell[N].hout;
       M_node_[1:N]=Cell[:].m_in;
       M_node_[N+1]=-Cell[N].m_out;//change sign
      algorithm
        Q:=0;
        for i in 1:N loop
         Q:=Q + Cell[i].Q;
        end for;
       annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
             Rectangle(
               extent={{-90,42},{92,-28}},
               lineColor={28,108,200},
               fillColor={170,255,255},
               fillPattern=FillPattern.Solid),
             Ellipse(
               extent={{-100,42},{-66,-28}},
               lineColor={28,108,200},
               fillColor={0,128,255},
               fillPattern=FillPattern.Solid),
             Ellipse(
               extent={{66,42},{100,-28}},
               lineColor={28,108,200},
               fillColor={0,128,255},
               fillPattern=FillPattern.Solid),
             Rectangle(
               extent={{-72,16},{40,-2}},
               lineColor={28,108,200},
               fillColor={0,127,0},
               fillPattern=FillPattern.Solid),
             Polygon(
               points={{20,26},{46,6},{22,-18},{66,6},{20,26}},
               lineColor={28,108,200},
               fillColor={0,127,0},
               fillPattern=FillPattern.Solid),
             Text(
               extent={{-46,40},{24,20}},
               lineColor={28,108,200},
               fillColor={0,127,0},
               fillPattern=FillPattern.Solid,
                textString="VeryFlow")}),                             Diagram(
             coordinateSystem(preserveAspectRatio=false)));
      end veryFlow;

      model CrossFlow_liquid

        TJUthermo.Interface.SimpleFlangeMulti node_a(N=N)
          annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
        TJUthermo.Interface.SimpleFlangeMulti node_b(N=N)
          annotation (Placement(transformation(extent={{88,-8},{108,12}})));
        Interface.mutiThermoP node_wall(N=N)
          annotation (Placement(transformation(extent={{-46,10},{40,52}})));
          //数据汇总//
        TJUthermo.DataRecord.Arrays Summary(n=N,T=Cell[:].T);
          //常量//
          parameter Integer N(min=2)=2;
          parameter Modelica.SIunits.Area F=0.0063;
          parameter Modelica.SIunits.Volume V=0.7263;
          Modelica.SIunits.MassFlowRate m_flow=0.5;
          parameter Modelica.SIunits.Temperature T1=313 annotation(Dialog(tab="Initialization"));
          parameter Modelica.SIunits.Temperature Tn=313 annotation(Dialog(tab="Initialization"));
          parameter Modelica.SIunits.Temperature Tstart[N]=linspace(T1,Tn,N) annotation(Dialog(tab="Initialization"));
          parameter Modelica.SIunits.CoefficientOfHeatTransfer  U;
          parameter Modelica.SIunits.SpecificHeatCapacity Cp;
          parameter Modelica.SIunits.Mass M;

          //变量//
          Modelica.SIunits.HeatFlowRate Q;
          replaceable TJUthermo.Components.units.Cell_liquid Cell[N](
          each F=F/N,
          each V=V/N,
          each U=U,
          each Cp=Cp,
          each M=M/N,
          Tstart=Tstart);

          TJUthermo.Interface.ThermoConvet ThermoC(final N=N);
          TJUthermo.Interface.nodeConvet_liquid Convet1(final N=N);
          TJUthermo.Interface.nodeConvet_liquid Convet2(final N=N);
      equation
        connect(node_a,Convet1.multi);
        connect(Convet1.single,Cell.A);
        connect(Cell.B,Convet2.single);
        connect(Convet2.multi,node_b);
        //连接壁面温度//
        connect(Cell.wall,ThermoC.single);
        connect(ThermoC.multi,node_wall);
      algorithm
        Q:=0;
        for i in 1:N loop
        Q:=Q + Cell[i].Q;
        end for;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-100,46},{100,-42}},
                lineColor={28,108,200},
                fillColor={170,170,255},
                fillPattern=FillPattern.Forward), Text(
                extent={{-68,10},{62,-30}},
                lineColor={28,108,200},
                fillColor={170,170,255},
                fillPattern=FillPattern.Solid,
                textString="LiquidFlow"),
              Rectangle(
                extent={{-86,18},{62,10}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={175,175,175}),
              Polygon(
                points={{58,28},{84,16},{60,2},{60,12},{58,28}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={175,175,175})}), Diagram(coordinateSystem(
                preserveAspectRatio=false), graphics={Text(
                extent={{-98,16},{-54,-16}},
                lineColor={28,108,200},
                textString="In"), Text(
                extent={{50,18},{94,-14}},
                lineColor={28,108,200},
                textString="out")}));
      end CrossFlow_liquid;

      model CrossLiquid
        parameter Integer N=4;
        TJUthermo.Interface.mutiThermoP wall(N=N)
          annotation (Placement(transformation(extent={{-14,8},{12,34}})));
          //常量区//
          constant Real pi=3.1415;
          parameter Modelica.SIunits.Area F=0.0063;
          parameter Modelica.SIunits.Volume V=0.7263;
          Modelica.SIunits.MassFlowRate m_flow[N];//修改处
          parameter Modelica.SIunits.Temperature Tstart=313 annotation(Dialog(tab="Initialization"));
          parameter Modelica.SIunits.CoefficientOfHeatTransfer  U;
          parameter Modelica.SIunits.SpecificHeatCapacity Cp;
          parameter Modelica.SIunits.Mass M;
          Modelica.SIunits.HeatFlowRate Q[N];//【w】
          Modelica.SIunits.Temperature Twall[N] "K";
          Modelica.SIunits.Temperature T[N];//(start=Tstart);
          Modelica.SIunits.Temperature Tin[N];//(start=Tstart);
          Modelica.SIunits.Temperature Tout[N];//(start=Tstart);
        Interface.SimpleFlangeMulti A(N=N)
          annotation (Placement(transformation(extent={{-108,-10},{-88,10}})));
        Interface.SimpleFlangeMulti B(N=N)
          annotation (Placement(transformation(extent={{88,-10},{108,10}})));
      initial equation
        for i in 1:N loop
        T[i]=Tstart;
        end for;
      equation
        for i in 1:N loop
        Twall[i]=wall.T[i];
        wall.Q[i]=Q[i];
        Q[i]=U*F/N*(Twall[i]-T[i]);
        m_flow[i]=A.m_flow[i];
        A.m_flow[i]+B.m_flow[i]=0;
        Tin[i]=A.T[i];
        B.T[i]=Tout[i];
        m_flow[i]*Cp*(Tin[i]-Tout[i])+Q[i]=M/N*Cp*der(T[i]);//2*h=hin+hout;
        Tout[i]=2*T[i]-Tin[i];
        end for;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-96,38},{96,-50}},
                lineColor={28,108,200},
                fillColor={85,170,255},
                fillPattern=FillPattern.Forward), Text(
                extent={{-58,6},{50,-38}},
                lineColor={28,108,200},
                fillColor={85,170,255},
                fillPattern=FillPattern.Forward,
                textString="CrossLiquid")}),                           Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end CrossLiquid;

      model Cell_collector_new

        TJUthermo.Interface.SimpleFlange_a inFlow
          annotation (Placement(transformation(extent={{-108,-6},{-88,14}})));
        TJUthermo.Interface.SimpleFlange_b outFlow
          annotation (Placement(transformation(extent={{88,-6},{108,14}})));
          parameter Modelica.SIunits.Area A;
          parameter Real K;
          parameter Real K_cp=0 "导热油比热斜率";
          parameter Real b_cp=2.2e3 "导热油比热截距";
          parameter Real K_rou=0 "导热油比热斜率";
          parameter Real b_rou=900 "导热油比热截距";
          parameter Modelica.SIunits.Mass M;
          Modelica.SIunits.SpecificHeatCapacity Cp;
          Modelica.SIunits.MassFlowRate m_flow;
          parameter Modelica.SIunits.Temperature T0 "初始温度";
          parameter Boolean use_eff2=true;
          parameter Modelica.SIunits.Conversions.NonSIunits.Time_hour ti=10
          "开始模拟太阳时刻";
          parameter Modelica.SIunits.Conversions.NonSIunits.Angle_deg lati=30
          "所在地纬度";
          parameter Integer Date=100 "日期";
          Modelica.SIunits.Temperature dt;
          Modelica.SIunits.Temperature Tin;
          Modelica.SIunits.Temperature Tout;
          Modelica.SIunits.Temperature T(start=T0);
          Modelica.SIunits.Temperature Tam;
          Modelica.SIunits.HeatFlowRate Q;//【w】
          Real eta;
        Modelica.SIunits.HeatFlux I;
                 //输入变量数//
        Modelica.Blocks.Interfaces.RealInput windspeed annotation (Placement(transformation(
              extent={{-17,-17},{17,17}},
              rotation=270,
              origin={-65,33}), iconTransformation(
              extent={{-9,-9},{9,9}},
              rotation=270,
              origin={-57,25})));
        Modelica.Blocks.Interfaces.RealInput I_in annotation (Placement(transformation(
              extent={{-20,-20},{20,20}},
              rotation=-90,
              origin={46,36}), iconTransformation(
              extent={{-9,-9},{9,9}},
              rotation=-90,
              origin={57,25})));
        Modelica.Blocks.Interfaces.RealInput Tam_in annotation (Placement(transformation(
                extent={{-28,22},{12,62}}), iconTransformation(
              extent={{-8,-8},{8,8}},
              rotation=-90,
              origin={-2,28})));
      //initial equation
        //T=T0;
      equation
          if cardinality(windspeed) == 0 then
          windspeed = 3 "风速";
          end if;
          if cardinality(Tam_in) == 0 then
          Tam_in = 25+273.15 "外界温度";
          end if;
          if cardinality(I_in) == 0 then
          I_in = 700 "外界辐射量";
          end if;
          if use_eff2 then
        eta=TJUthermo.HeatTransferModel.solar_PTC(dt=dt,ti=time,t0=10,lati=lati,N=Date,I=I);
          else
        eta=(K*(73.6-0.004206*dt)+7.44*(dt/I)-0.00958*(dt*dt/I))*0.01;
          end if;
       Cp=1000*(0.0022*T+1.3192);
       Tin=inFlow.T;
       Tout=outFlow.T;
       Tam=Tam_in+273.15;
       dt=0.5*(Tin+Tout)-Tam;
       I=I_in;
       Q=K*A*I*eta;
       Q+m_flow*Cp*(Tin-Tout)=M*Cp*der(T);
       m_flow=inFlow.m_flow;
       inFlow.m_flow+outFlow.m_flow=0;
       Tin+Tout=2*T;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-98,28},{98,-22}},
                lineColor={28,108,200},
                fillColor={255,255,170},
                fillPattern=FillPattern.Solid), Text(
                extent={{-80,18},{78,-12}},
                lineColor={28,108,200},
                fillColor={255,255,170},
                fillPattern=FillPattern.Solid,
                textString="collector unit"),
              Text(
                extent={{-74,22},{-60,18}},
                lineColor={28,108,200},
                fillColor={255,255,170},
                fillPattern=FillPattern.Solid,
                textString="风速"),
              Text(
                extent={{42,22},{52,16}},
                lineColor={28,108,200},
                fillColor={255,255,170},
                fillPattern=FillPattern.Solid,
                textString="辐射量"),
              Text(
                extent={{8,24},{26,16}},
                lineColor={28,108,200},
                fillColor={255,255,170},
                fillPattern=FillPattern.Solid,
                textString="外界温度"),
              Text(
                extent={{-36,-12},{38,-40}},
                lineColor={28,108,200},
                textString="NEW")}),            Diagram(coordinateSystem(
                preserveAspectRatio=false)));
      end Cell_collector_new;
    end units;

    package main_components
      model Hx1counter_liquid

         parameter Integer N(min=1)=2 "分割单元体数目";
         parameter Modelica.SIunits.CoefficientOfHeatTransfer U1=800
          "流体换热系数";
         parameter Modelica.SIunits.CoefficientOfHeatTransfer U2=800
          "流体换热系数";
         parameter Modelica.SIunits.Area F_liquid1=6.3 "工质换热面积";
         parameter Modelica.SIunits.Area F_liquid2=6.3 "流体换热面积";
         parameter Modelica.SIunits.Volume V_liquid1=0.7263 "流体容积";
         parameter Modelica.SIunits.Volume V_liquid2=0.7263 "流体容积";
         parameter Modelica.SIunits.SpecificHeatCapacity Cp1 "流体比热";
         parameter Modelica.SIunits.SpecificHeatCapacity Cp2 "流体比热";
         parameter Modelica.SIunits.Mass M1 "流体质量";
         parameter Modelica.SIunits.Mass M2 "流体质量";
         parameter Modelica.SIunits.Mass M_wall "管壁质量";
         parameter Modelica.SIunits.SpecificHeatCapacity Cp_wall "管壁比热";
         Modelica.SIunits.HeatFlowRate Q;
         //==============================初值条件区===================================//
         parameter Modelica.SIunits.Temperature T1_wall=0.5*(T2n_liquid+T1_liquid)
          "墙壁初始温度"                                                                            annotation(Dialog(tab="初始化参数"));
         parameter Modelica.SIunits.Temperature Tn_wall=0.5*(T1n_liquid+T2_liquid)
          "墙壁初始温度"                                                                             annotation(Dialog(tab="初始化参数"));
         parameter Modelica.SIunits.Temperature T1_liquid "流体入口温度" annotation(Dialog(tab="初始化参数"));
         parameter Modelica.SIunits.Temperature T1n_liquid "流体出口温度"
                                                                                annotation(Dialog(tab="初始化参数"));
         parameter Modelica.SIunits.Temperature T2_liquid "流体入口温度" annotation(Dialog(tab="初始化参数"));
         parameter Modelica.SIunits.Temperature T2n_liquid "流体出口温度"
                                                                                annotation(Dialog(tab="初始化参数"));
         //=========================================================================//
        units.SimpleFlowliquid_type2 simpleFlowliquid1(
          final N=N,
          final F=F_liquid1,
          final V=V_liquid1,
          final T1=T1_liquid,
          final Tn=T1n_liquid,
          final U=U1,
          final Cp=Cp1,
          final M=M1)
          annotation (Placement(transformation(extent={{-58,-60},{54,-22}})));
        units.SimpleFlowliquid_type2 simpleFlowliquid2(
          final N=N,
          final F=F_liquid2,
          final V=V_liquid2,
          final T1=T2_liquid,
          final Tn=T2n_liquid,
          final U=U2,
          final Cp=Cp2,
          final M=M2)
          annotation (Placement(transformation(extent={{52,68},{-58,30}})));
        source_sink_ports.CounterDoubleWalls Walls(N=N,
          final T1=T1_wall,
          final Tn=Tn_wall,
          final M=M_wall,
          final Cp=Cp_wall)
          annotation (Placement(transformation(extent={{-36,-8},{26,26}})));
        Interface.SimpleFlange_a inflow1
          annotation (Placement(transformation(extent={{-112,-50},{-92,-30}})));
        Interface.SimpleFlange_a inflow2
          annotation (Placement(transformation(extent={{88,38},{108,58}})));
        Interface.SimpleFlange_b outflow2
          annotation (Placement(transformation(extent={{-110,40},{-90,60}})));
        Interface.SimpleFlange_b outflow1 annotation (Placement(transformation(extent=
                 {{88,-50},{108,-30}}), iconTransformation(extent={{88,-50},{108,-30}})));
      equation
        Q=simpleFlowliquid1.Q;
        connect(simpleFlowliquid2.node_wall, Walls.Wall_side1) annotation (Line(
              points={{-1.35,43.11},{-1.35,28.555},{-5,28.555},{-5,13.42}}, color={28,
                108,200}));
        connect(Walls.Wall_side2, simpleFlowliquid1.node_wall) annotation (Line(
              points={{-4.38,2.37},{-4.38,-16.815},{-3.68,-16.815},{-3.68,-35.11}},
              color={28,108,200}));
        connect(inflow1, simpleFlowliquid1.node_a) annotation (Line(points={{-102,-40},
                {-58,-40},{-58,-41}}, color={28,108,200}));
        connect(simpleFlowliquid2.node_a, inflow2) annotation (Line(points={{52,49},{74,
                49},{74,48},{98,48}}, color={28,108,200}));
        connect(simpleFlowliquid2.node_b, outflow2) annotation (Line(points={{-56.9,48.62},
                {-56.9,49.31},{-100,49.31},{-100,50}}, color={28,108,200}));
        connect(simpleFlowliquid1.node_b, outflow1) annotation (Line(points={{52.88,-40.62},
                {77.44,-40.62},{77.44,-40},{98,-40}}, color={28,108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
                extent={{-100,78},{100,-76}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-100,66},{100,26}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Forward),
              Rectangle(
                extent={{-100,-22},{100,-62}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Forward),
              Polygon(
                points={{50,-56},{76,-40},{50,-28},{50,-56}},
                lineColor={28,108,200},
                fillColor={170,213,255},
                fillPattern=FillPattern.Forward),
              Rectangle(
                extent={{-82,-38},{50,-46}},
                lineColor={28,108,200},
                fillColor={170,213,255},
                fillPattern=FillPattern.Forward),
              Rectangle(
                extent={{-54,50},{82,42}},
                lineColor={28,108,200},
                fillColor={170,213,255},
                fillPattern=FillPattern.Forward),
              Polygon(
                points={{-52,32},{-78,48},{-52,60},{-52,32}},
                lineColor={28,108,200},
                fillColor={170,213,255},
                fillPattern=FillPattern.Forward),
              Text(
                extent={{-80,24},{80,-18}},
                lineColor={28,108,200},
                fillColor={170,213,255},
                fillPattern=FillPattern.Forward,
                textString="Hx1counter-liquid"),
              Text(
                extent={{-30,-30},{8,-58}},
                lineColor={28,108,200},
                fillColor={170,213,255},
                fillPattern=FillPattern.Forward,
                textString="1"),
              Text(
                extent={{-22,60},{10,32}},
                lineColor={28,108,200},
                fillColor={170,213,255},
                fillPattern=FillPattern.Forward,
                textString="2")}),                 Diagram(coordinateSystem(
                preserveAspectRatio=false)));
      end Hx1counter_liquid;

      model Hx1simple
      replaceable package Medium =TJUthermo.Media constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model" annotation(choicesAllMatching=true);
        Interface.Flange_a InFlow_working
          annotation (Placement(transformation(extent={{-110,-40},{-90,-20}})));
        Interface.Flange_b OutFlow_working
          annotation (Placement(transformation(extent={{88,-40},{108,-20}})));
        Interface.SimpleFlange_a InFlow_liquid
          annotation (Placement(transformation(extent={{-110,28},{-90,48}})));
        Interface.SimpleFlange_b OutFlow_liquid
          annotation (Placement(transformation(extent={{88,28},{108,48}})));
          //===============常量================//
           parameter Integer N(min=1)=2 "分割单元体数目";
           parameter Modelica.SIunits.CoefficientOfHeatTransfer Ul=1000
          "液态换热系数";
            parameter Modelica.SIunits.CoefficientOfHeatTransfer Utp=5000
          "两相换热系数";
            parameter Modelica.SIunits.CoefficientOfHeatTransfer Uv=2000
          "气态换热系数";
            parameter Modelica.SIunits.CoefficientOfHeatTransfer U=800
          "流体换热系数";
            parameter Modelica.SIunits.Area F_working=6.3 "工质换热面积";
            parameter Modelica.SIunits.Area F_liquid=6.3 "流体换热面积";
            parameter Modelica.SIunits.Volume V_working=0.7263 "工质容积";
            parameter Modelica.SIunits.Volume V_liquid=0.7263 "流体容积";
            parameter Modelica.SIunits.Mass M "流体质量";
            parameter Modelica.SIunits.SpecificHeatCapacity Cp "流体比热";
            parameter Modelica.SIunits.Mass M_wall "管壁质量";
            parameter Modelica.SIunits.SpecificHeatCapacity Cp_wall
          "管壁比热";
            parameter Real max_drhode=100
          "控制密度随压力焓值变化的偏导数值，推荐取100";
            parameter Real witdth_x=0.05 "减缓突变设定值";
            //===============初值区=================//
            parameter Modelica.SIunits.Temperature Tin_start=303 annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.Temperature Tout_start annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.Temperature T1_wall=0.5*(Tin_start+T1_liquid)
          "墙壁初始温度"                                                                            annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.Temperature Tn_wall=0.5*(Tout_start+Tn_liquid)
          "墙壁初始温度"                                                                             annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.Temperature T1_liquid
          "流体入口温度"                                                    annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.Temperature Tn_liquid
          "流体出口温度"                                                    annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.MassFlowRate m_working=0.5
          "工质质量流量"                                                         annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.MassFlowRate m_liquid=0.5
          "流体质量流量"                                                        annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.Pressure pstart=1.2e6 annotation(Dialog(tab="初始化参数"));
          //===============变量================//
          replaceable TJUthermo.Components.units.SimpleFlow Flow1(
          redeclare final package Medium=Medium,
          final N=N,
          final Ul=Ul,
          final Utp=Utp,
          final Uv=Uv,
          final F=F_working,
          final V=V_working,
          final max_drdp=max_drhode,
          final max_drdh=max_drhode,
          final witdth_x=witdth_x,
          final Tin_start=Tin_start,
          final Tout_start=Tout_start,
          final m_flow_start=m_working,
          final pstart=pstart);
          replaceable TJUthermo.Components.units.SimpleFlowliquid FlowLiquids(
          final N=N,
          final F=F_liquid,
          final V=V_liquid,
          final m_flow=m_liquid,
          final T1=T1_liquid,
          final Tn=Tn_liquid,
          final U=U,
          final Cp=Cp,
          final M=M);

          TJUthermo.Components.source_sink_ports.DoubleWalls CWalls(N=N,
          final T1=T1_wall,
          final Tn=Tn_wall,
          final M=M_wall,
          final Cp=Cp_wall);
      equation
          connect(InFlow_working,Flow1.InFlow);
          connect(Flow1.OutFlow,OutFlow_working);
          connect(InFlow_liquid,FlowLiquids.node_a);
          connect(FlowLiquids.node_b,OutFlow_liquid);
          connect(Flow1.Wall,CWalls.Wall_side1);
          connect(FlowLiquids.node_wall,CWalls.Wall_side2);
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
                extent={{-100,64},{100,-66}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={175,175,175}),
              Text(
                extent={{-82,22},{72,-18}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={128,0,255},
                textString="Hx1—Simple"),
              Rectangle(
                extent={{-88,44},{86,28}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={128,0,255}),
              Rectangle(
                extent={{-88,-24},{86,-40}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={128,0,255}),
              Text(
                extent={{-62,-70},{56,-98}},
                lineColor={28,108,200},
                textString="顺流")}),                                Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
                extent={{-88,46},{86,30}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={128,0,255}),
              Rectangle(
                extent={{-88,-22},{86,-38}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={128,0,255}),
              Text(
                extent={{-82,24},{72,-16}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={128,0,255},
                textString="Hx1—Simple")}));
      end Hx1simple;

      model Hx1counter
        //========================工质=======================//
        replaceable package Medium =TJUthermo.Media constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model" annotation(choicesAllMatching=true);
      //=========================接口区============================//
        Interface.Flange_a InFlow_working
          annotation (Placement(transformation(extent={{-110,-44},{-90,-24}})));
        Interface.Flange_b OutFlow_working
          annotation (Placement(transformation(extent={{88,-42},{108,-22}})));
        Interface.SimpleFlange_a InFlow_liquid
          annotation (Placement(transformation(extent={{88,36},{108,56}})));
        Interface.SimpleFlange_b OutFlow_liquid
          annotation (Placement(transformation(extent={{-110,34},{-90,54}})));
           //===============常量================//
           parameter Integer N(min=1)=2 "分割单元体数目";
           parameter Modelica.SIunits.CoefficientOfHeatTransfer Ul=1000
          "液态换热系数";
            parameter Modelica.SIunits.CoefficientOfHeatTransfer Utp=5000
          "两相换热系数";
            parameter Modelica.SIunits.CoefficientOfHeatTransfer Uv=2000
          "气态换热系数";
            parameter Modelica.SIunits.CoefficientOfHeatTransfer U=800
          "流体换热系数";
            parameter Boolean Mcons=false;
            parameter Modelica.SIunits.Area F_working=6.3 "工质换热面积";
            parameter Modelica.SIunits.Area F_liquid=6.3 "流体换热面积";
            parameter Modelica.SIunits.Volume V_working=0.7263 "工质容积";
            parameter Modelica.SIunits.Volume V_liquid=0.7263 "流体容积";
            parameter Modelica.SIunits.Mass M "流体质量";
            parameter Modelica.SIunits.SpecificHeatCapacity Cp "流体比热";
            parameter Modelica.SIunits.Mass M_wall "管壁质量";
            parameter Modelica.SIunits.SpecificHeatCapacity Cp_wall
          "管壁比热";
            parameter Real max_drhode=100
          "控制密度随压力焓值变化的偏导数值，推荐取100";
            parameter Real witdth_x=0.05 "减缓突变设定值";
            //===============初值区=================//
            parameter Modelica.SIunits.Temperature Tin_start=303 annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.Temperature Tout_start annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.Temperature T1_wall=0.5*(Tout_start+T1_liquid)
          "墙壁初始温度"                                                                            annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.Temperature Tn_wall=0.5*(Tin_start+Tn_liquid)
          "墙壁初始温度"                                                                             annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.Temperature T1_liquid
          "流体入口温度"                                                    annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.Temperature Tn_liquid
          "流体出口温度"                                                    annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.MassFlowRate m_working=0.5
          "工质质量流量"                                                         annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.MassFlowRate m_liquid=0.5
          "流体质量流量"                                                        annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.Pressure pstart=1.2e6 annotation(Dialog(tab="初始化参数"));
            Modelica.SIunits.HeatFlowRate Q;
          //===============变量================//
        source_sink_ports.CounterDoubleWalls Walls(N=N,
          final T1=T1_wall,
          final Tn=Tn_wall,
          final M=M_wall,
          final Cp=Cp_wall)
          annotation (Placement(transformation(extent={{-16,-8},{4,12}})));
        units.SimpleFlow workingFlow(
          redeclare final package Medium=Medium,
          final N=N,
          final Ul=Ul,
          final Utp=Utp,
          final Uv=Uv,
          final Mcons=Mcons,
          final F=F_working,
          final V=V_working,
          final max_drdp=max_drhode,
          final max_drdh=max_drhode,
          final witdth_x=witdth_x,
          final Tin_start=Tin_start,
          final Tout_start=Tout_start,
          final m_flow_start=m_working,
          final pstart=pstart)
          annotation (Placement(transformation(extent={{-32,-52},{38,-16}})));
        units.SimpleFlowliquid Flowliquid(
          final N=N,
          final F=F_liquid,
          final V=V_liquid,
          final T1=T1_liquid,
          final Tn=Tn_liquid,
          final U=U,
          final Cp=Cp,
          final M=M)
          annotation (Placement(transformation(extent={{38,64},{-34,26}})));

      equation
        Q=workingFlow.Q;
        connect(Flowliquid.node_wall, Walls.Wall_side1) annotation (Line(points={{3.08,
                39.11},{3.08,21.555},{-6,21.555},{-6,4.6}}, color={28,108,200}));
        connect(Walls.Wall_side2, workingFlow.Wall) annotation (Line(points={{-5.8,-1.9},
                {-5.8,-14.95},{3.7,-14.95},{3.7,-25.36}}, color={28,108,200}));
        connect(InFlow_working, workingFlow.InFlow) annotation (Line(points={{-100,-34},
                {-30.6,-34},{-30.6,-33.28}}, color={0,0,255}));
        connect(workingFlow.OutFlow, OutFlow_working) annotation (Line(points={{37.3,-32.92},
                {68.65,-32.92},{68.65,-32},{98,-32}}, color={0,0,255}));
        connect(OutFlow_liquid, Flowliquid.node_b) annotation (Line(points={{-100,44},
                {-33.28,44},{-33.28,44.62}}, color={28,108,200}));
        connect(Flowliquid.node_a, InFlow_liquid) annotation (Line(points={{38,45},{68,
                45},{68,46},{98,46}}, color={28,108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
                extent={{-100,68},{100,-62}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={255,255,0}),
              Rectangle(
                extent={{-86,52},{80,34}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={255,0,0}),
              Rectangle(
                extent={{-84,-24},{82,-42}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={255,0,0}),
              Text(
                extent={{-78,24},{76,-10}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={255,0,0},
                textString="CounterFlow")}),                           Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end Hx1counter;

      model Pump_simple
      replaceable package medium =TJUthermo.Media constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model"
        annotation(choicesAllMatching=true);
        Interface.Flange_a InFlow
          annotation (Placement(transformation(extent={{-52,-34},{-32,-14}})));
        Interface.Flange_b OutFlow
          annotation (Placement(transformation(extent={{22,36},{42,56}})));
          parameter Real n_pump(min=0,max=1)=0.67;
          parameter Modelica.SIunits.MassFlowRate m_flow=0.5;
          parameter Modelica.SIunits.Pressure p0=0.3e6;
          parameter Modelica.SIunits.Pressure pstart_out=1.2e6;
          parameter Modelica.SIunits.Temperature T0;
          parameter Boolean steadyInFlow=true
          "进口在一段时间内设为稳态";
          parameter Modelica.SIunits.Time steadytime=5
          "开始稳定段设置时间";
          Modelica.SIunits.SpecificEnthalpy hin(start=medium.specificEnthalpy_pT(p0,T0));
          Modelica.SIunits.SpecificEnthalpy hout(start=medium.specificEnthalpy_pT(pstart_out,T0));
          Modelica.SIunits.SpecificEnthalpy hout_s(start=medium.specificEnthalpy_pT(p0,T0));
          medium.ThermodynamicState state_in;
          medium.ThermodynamicState state_out;
          Modelica.SIunits.Pressure pin(start=p0);
          Modelica.SIunits.Pressure pout(start=pstart_out);
          Modelica.SIunits.SpecificEntropy s_in;
      //initial equation
        //hin=medium.specificEnthalpy_pT(p0,T0);
        //state_in=medium.setState_pT(p0,T0);
        //state_out=medium.setState_pT(pstart_out,T0);
      equation
        //InFlow.m_flow=m_flow;
        OutFlow.m_flow=-m_flow;
        if steadyInFlow then
          if time<steadytime then
          hin=medium.specificEnthalpy_pT(p0,T0);
          else
        hin=inStream(InFlow.h);
          end if;
        else
        hin=inStream(InFlow.h);
        end if;
        InFlow.h=hin;
        hout=OutFlow.h;
        //OutFlow.h=inStream(OutFlow.h);//修改//
        pin=InFlow.p;
        pin=p0;
        pout=OutFlow.p;
        state_in=medium.setState_ph(pin,hin);
        state_out=medium.setState_ph(pout,hout);
        s_in=state_in.s;
        hout_s=medium.specificEnthalpy_ps(pout,s_in);
        (hout_s-hin)/(hout-hin)=n_pump;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
                extent={{-44,48},{32,-28}},
                lineColor={28,108,200},
                fillColor={85,255,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{8,38},{-28,10},{8,-22},{8,-20},{8,38}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-18,18},{4,-2}},
                lineColor={28,108,200},
                fillColor={85,255,255},
                fillPattern=FillPattern.Solid,
                textString="Pump")}), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end Pump_simple;

      model Pump_type2
      replaceable package medium =TJUthermo.Media constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model"
        annotation(choicesAllMatching=true);
        Interface.Flange_a InFlow
          annotation (Placement(transformation(extent={{-52,-34},{-32,-14}})));
        Interface.Flange_b OutFlow
          annotation (Placement(transformation(extent={{22,36},{42,56}})));
          parameter Real n_pump(min=0,max=1)=0.67;
          Modelica.SIunits.MassFlowRate m_flow;
          parameter Modelica.SIunits.Pressure p0=0.3e6;
          parameter Modelica.SIunits.Pressure pstart_out=1.2e6;
          parameter Modelica.SIunits.Temperature T0;
          parameter Boolean steadyInFlow=false
          "进口在一段时间内设为稳态";
          parameter Modelica.SIunits.Time steadytime=5
          "开始稳定段设置时间";
          parameter Modelica.SIunits.Volume v_s;
          parameter Modelica.SIunits.Frequency Np=30;
          Modelica.SIunits.Density rou_in;
          Modelica.SIunits.SpecificEnthalpy hin(start=medium.specificEnthalpy_pT(p0,T0));
          Modelica.SIunits.SpecificEnthalpy hout(start=medium.specificEnthalpy_pT(pstart_out,T0));
          Modelica.SIunits.SpecificEnthalpy hout_s(start=medium.specificEnthalpy_pT(p0,T0));
          //medium.ThermodynamicState state_in;
          //medium.ThermodynamicState state_out;
          Modelica.SIunits.Pressure pin(start=p0);
          Modelica.SIunits.Pressure pout(start=pstart_out);
          Modelica.SIunits.SpecificEntropy s_in;
      //initial equation
        //hin=medium.specificEnthalpy_pT(p0,T0);
        //state_in=medium.setState_pT(p0,T0);
        //state_out=medium.setState_pT(pstart_out,T0);
      equation
        //InFlow.m_flow=m_flow;
        m_flow=v_s*rou_in*Np;
        OutFlow.m_flow=-m_flow;
        if steadyInFlow then
          if time<steadytime then
          hin=medium.specificEnthalpy_pT(p0,T0);
          else
        hin=inStream(InFlow.h);
          end if;
        else
        hin=inStream(InFlow.h);
        end if;
        InFlow.h=hin;
        hout=OutFlow.h;
        //OutFlow.h=inStream(OutFlow.h);//修改//
        pin=InFlow.p;
        m_flow=InFlow.m_flow;
        pout=OutFlow.p;
        //state_in=medium.setState_ph(pin,hin);
        rou_in=medium.density_ph(pin,hin);
        //state_out=medium.setState_ph(pout,hout);
        s_in=medium.specificEntropy_ph(pin,hin);
        hout_s=medium.specificEnthalpy_ps(pout,s_in);
        (hout_s-hin)/(hout-hin)=n_pump;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
                extent={{-44,48},{32,-28}},
                lineColor={28,108,200},
                fillColor={85,255,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-8,40},{-44,12},{-8,-20},{-8,-18},{-8,40}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-34,22},{-12,2}},
                lineColor={28,108,200},
                fillColor={85,255,255},
                fillPattern=FillPattern.Solid,
                textString="Pump"),
              Rectangle(
                extent={{-8,30},{50,-6}},
                lineColor={28,108,200},
                fillColor={0,128,255},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-2,24},{42,0}},
                lineColor={255,255,255},
                fillColor={0,128,255},
                fillPattern=FillPattern.Solid,
                textString="type 2")}),
                                      Diagram(coordinateSystem(preserveAspectRatio=false)));
      end Pump_type2;

      model Turbine_simple
      replaceable package medium =TJUthermo.Media constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model"
        annotation(choicesAllMatching=true);
        Interface.Flange_a InFlow
          annotation (Placement(transformation(extent={{-60,32},{-40,52}})));
        Interface.Flange_b OutFlow
          annotation (Placement(transformation(extent={{36,-76},{56,-56}})));
          parameter Real n_turbine(min=0,max=1)=0.67;
          parameter Modelica.SIunits.Volume V_s=0.5;
          Modelica.SIunits.Density rou_in;
          parameter Modelica.SIunits.Frequency Np=30;
          parameter Real n_ele=0.88;
          parameter Modelica.SIunits.Pressure pin_start=1.2e6 annotation(Dialog(tab="Initialization"));
          parameter Modelica.SIunits.Pressure pout_start=0.3e6 annotation(Dialog(tab="Initialization"));
          parameter Modelica.SIunits.Temperature Tin_start=110+273 annotation(Dialog(tab="Initialization"));
         // parameter medium.ThermodynamicState state0=medium.setState_pT(pin_start,Tin_start);
          parameter Modelica.SIunits.SpecificEntropy s_start=medium.specificEntropy_pT(pin_start,Tin_start)
                                                                                                    annotation(Dialog(tab="Initialization"));
          Modelica.SIunits.MassFlowRate m_flow;
          Modelica.SIunits.SpecificEnthalpy hin;//(start=medium.specificEnthalpy_pT(pin_start,Tin_start));
          Modelica.SIunits.SpecificEnthalpy hout;
          Modelica.SIunits.SpecificEnthalpy hout_s;
          medium.ThermodynamicState state_in(phase(
                                             start =     1),s(start=s_start));
                                                            //(start=state0);
          medium.ThermodynamicState state_out;
          Modelica.SIunits.Pressure pin(start=pin_start);
          Modelica.SIunits.Pressure pout;//(start=pout_start);
          Modelica.SIunits.SpecificEntropy s_in(start=s_start);
          Modelica.SIunits.Power W_turbine;
      //initial equation
        //pin=pin_start;
        //hin=medium.specificEnthalpy_pT(pin_start,Tin_start);
          //state_in=medium.setState_pT(pin_start,Tin_start);
          //state_out=medium.setState_ps(pout_start,s_start);
          //hout=medium.specificEnthalpy_ps(pout_start,s_start);
         //hout_s=medium.specificEnthalpy_ps(pout_start,s_start);
      equation
        m_flow=V_s*rou_in*Np;
        InFlow.m_flow=m_flow;
        OutFlow.m_flow=-m_flow;
        hin=inStream(InFlow.h);
        hout=OutFlow.h;
        InFlow.h=hin;
        pin=InFlow.p;
        pout=OutFlow.p;
        state_in=medium.setState_ph(pin,hin);
        rou_in=state_in.d;
        state_out=medium.setState_ph(pout,hout);
        s_in=state_in.s;
        hout_s=medium.specificEnthalpy_ps(pout,s_in);
        (hout-hin)/(hout_s-hin)=n_turbine;
        W_turbine=m_flow*(hin-hout)*n_ele;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Polygon(
                points={{-52,42},{46,88},{46,-68},{-52,-22},{-52,42}},
                lineColor={28,108,200},
                fillColor={85,255,255},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-40,34},{36,26}},
                lineColor={28,108,200},
                fillColor={0,255,0},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-38,30},{32,-14}},
                lineColor={28,108,200},
                fillColor={85,255,255},
                fillPattern=FillPattern.Solid,
                textString="Turbine"),
              Rectangle(
                extent={{-40,-6},{36,-14}},
                lineColor={28,108,200},
                fillColor={0,255,0},
                fillPattern=FillPattern.Solid)}),                      Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={
              Polygon(
                points={{-50,40},{48,86},{48,-70},{-50,-24},{-50,40}},
                lineColor={28,108,200},
                fillColor={85,255,255},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-36,30},{34,-14}},
                lineColor={28,108,200},
                fillColor={85,255,255},
                fillPattern=FillPattern.Solid,
                textString="Turbine"),
              Rectangle(
                extent={{-38,34},{38,26}},
                lineColor={28,108,200},
                fillColor={0,255,0},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-38,-6},{38,-14}},
                lineColor={28,108,200},
                fillColor={0,255,0},
                fillPattern=FillPattern.Solid)}));
      end Turbine_simple;

      model Tank
      replaceable package medium =TJUthermo.Media constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model"
        annotation(choicesAllMatching=true);
        Interface.Flange_a InFlow
          annotation (Placement(transformation(extent={{50,30},{70,50}})));
        Interface.Flange_b OutFlow
          annotation (Placement(transformation(extent={{-64,-60},{-44,-40}})));
          Modelica.SIunits.MassFlowRate min;
          Modelica.SIunits.MassFlowRate mout;
          parameter Modelica.SIunits.Pressure p0=0.3e6;
          parameter Modelica.SIunits.Temperature T0=40+273.15;
          parameter Modelica.SIunits.Area Aera=0.5;
          parameter Modelica.SIunits.Length L0=0.5;
          Modelica.SIunits.Length L;
          Modelica.SIunits.Density rou(start=medium.density_pT(p0,T0));
          Modelica.SIunits.SpecificEnthalpy h;
          Modelica.SIunits.SpecificEnthalpy hin;
          medium.ThermodynamicState State_tank;
      initial equation
        L=L0;
        h=medium.specificEnthalpy_pT(p0,T0);
      equation
        min=InFlow.m_flow;
        mout=OutFlow.m_flow;
        InFlow.p=p0;
        OutFlow.p=p0;
        hin=inStream(InFlow.h);
        hin=InFlow.h;
        OutFlow.h=h;
          min+mout=Aera*rou*der(L);
          min*hin+mout*h=Aera*rou*der(h*L);
          rou=medium.density_ph(p0,h);//恒压
          State_tank=medium.setState_ph(p0,h);

        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-52,66},{60,-72}},
                lineColor={28,108,200},
                fillColor={170,170,255},
                fillPattern=FillPattern.VerticalCylinder), Text(
                extent={{-36,30},{40,-36}},
                lineColor={28,108,200},
                fillPattern=FillPattern.VerticalCylinder,
                fillColor={170,170,255},
                textString="Tank")}),                                  Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end Tank;

      model SolarCollector

        Interface.SimpleFlange_a inFlow
          annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
        Interface.SimpleFlange_b outFlow
          annotation (Placement(transformation(extent={{90,-8},{110,12}})));
        Modelica.Blocks.Interfaces.RealInput windspeed annotation (Placement(transformation(
              extent={{-20,-20},{20,20}},
              rotation=-90,
              origin={-78,32}), iconTransformation(
              extent={{-8,-8},{8,8}},
              rotation=-90,
              origin={-72,58})));
        Modelica.Blocks.Interfaces.RealInput I_in annotation (Placement(transformation(
              extent={{-20,-20},{20,20}},
              rotation=-90,
              origin={-18,34}), iconTransformation(
              extent={{-8,-8},{8,8}},
              rotation=-90,
              origin={-12,58})));
        Modelica.Blocks.Interfaces.RealInput Tam_in annotation (Placement(transformation(
              extent={{-20,-20},{20,20}},
              rotation=-90,
              origin={34,70}), iconTransformation(
              extent={{-8,-8},{8,8}},
              rotation=-90,
              origin={46,58})));
          TJUthermo.DataRecord.Arrays Summary(n=N,T=Cell[:].T);
          Modelica.SIunits.HeatFlowRate Q;
          parameter Integer N(min=2) "划分单元体网格数";
          parameter Modelica.SIunits.Area Area "集热面积";
          parameter Real K "吸收系数";
          parameter Modelica.SIunits.Mass M "槽子中导热油质量";
          parameter Modelica.SIunits.SpecificHeatCapacity Cp "导热油比热";
          parameter Modelica.SIunits.MassFlowRate m_flow=1
          "导热油质量流量";
          parameter Modelica.SIunits.Temperature T0 "初始温度";
        TJUthermo.Components.units.Cell_collector Cell[N](
        each A=Area/N,
        each K=K,
        each M=M/N,
        each Cp=Cp,
        each T0=T0);

      equation
        if cardinality(windspeed) == 0 then
          windspeed=2.5;
        end if;
         connect(inFlow,Cell[1].inFlow);
         for i in 1:N-1 loop
           connect(Cell[i].outFlow,Cell[i+1].inFlow);
         end for;
         connect(Cell[N].outFlow,outFlow);
         for i in 1:N loop
         connect(I_in,Cell[i].I_in);
         connect(Tam_in,Cell[i].Tam_in);
         //connect(windspeed,Cell[i].windspeed);
         end for;
      algorithm
        Q:=0;
        for i in 1:N loop
        Q:=Q + Cell[i].Q;
        end for;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Polygon(
                points={{-96,58},{-86,48},{-78,32},{-76,20},{-72,-6},{-76,-36},{-76,-38},
                    {-86,-50},{-94,-56},{76,-58},{86,-48},{92,-30},{96,-16},{96,0},{96,
                    16},{90,38},{82,52},{78,58},{-96,58}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-100,18},{100,-16}},
                lineColor={28,108,200},
                fillColor={255,255,0},
                fillPattern=FillPattern.HorizontalCylinder),
              Text(
                extent={{-62,46},{46,24}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={255,255,0},
                textString="槽式"),
              Text(
                extent={{-56,-20},{52,-48}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={255,255,0},
                textString="集热器"),
              Text(
                extent={{-76,48},{-54,40}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={255,255,0},
                textString="风速"),
              Text(
                extent={{-24,80},{6,68}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={255,255,0},
                textString="辐射"),
              Text(
                extent={{54,54},{78,40}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={255,255,0},
                textString="温度")}),    Diagram(coordinateSystem(
                preserveAspectRatio=false), graphics={Rectangle(
                extent={{-100,26},{100,-24}},
                lineColor={28,108,200},
                fillColor={255,255,0},
                fillPattern=FillPattern.Solid)}));
      end SolarCollector;

      model Valve

        Interface.SimpleFlange_a inFlow
          annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
        Interface.SimpleFlange_b outFlow_1
          annotation (Placement(transformation(extent={{90,-8},{110,12}})));
        Interface.SimpleFlange_b outFlow_2 annotation (Placement(transformation(
                extent={{-10,-80},{10,-60}}), iconTransformation(extent={{-10,-80},{10,
                  -60}})));

       Real ratio(min=0,max=1);
        Modelica.Blocks.Interfaces.RealInput u annotation (Placement(transformation(
                extent={{-32,16},{8,56}}), iconTransformation(
              extent={{-20,-20},{20,20}},
              rotation=-90,
              origin={0,56})));
      equation
        ratio=u;
        inFlow.m_flow+outFlow_1.m_flow+outFlow_2.m_flow=0;
        ratio*inFlow.m_flow=-outFlow_1.m_flow;
        inFlow.T=outFlow_1.T;
        inFlow.T=outFlow_2.T;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
                extent={{-100,70},{100,-72}},
                lineColor={28,108,200},
                fillPattern=FillPattern.Solid,
                fillColor={215,215,215}),
              Polygon(
                points={{-86,54},{-86,-58},{-2,0},{-86,54}},
                lineColor={28,108,200},
                fillColor={0,255,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-2,0},{80,52},{80,-58},{-2,0}},
                lineColor={28,108,200},
                fillColor={0,255,255},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-26,-22},{18,-58}},
                lineColor={28,108,200},
                fillColor={0,255,255},
                fillPattern=FillPattern.Solid,
                textString="阀门")}), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end Valve;

      model Pump_liquid

        Interface.SimpleFlange_a inFlow
          annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
        Interface.SimpleFlange_b outFlow
          annotation (Placement(transformation(extent={{90,-10},{110,10}})));
         Modelica.SIunits.MassFlowRate m_flow;
         parameter Modelica.SIunits.MassFlowRate m0=1;
        Modelica.Blocks.Interfaces.RealInput u annotation (Placement(transformation(
              extent={{-20,-20},{20,20}},
              rotation=-90,
              origin={-2,26}), iconTransformation(
              extent={{-13,-13},{13,13}},
              rotation=-90,
              origin={1,39})));
      equation
      if cardinality(u) == 0 then
          m_flow = m0 "流速";
      else
          m_flow=u;
          end if;
        //inFlow.m_flow=m_flow;
        outFlow.m_flow=-m_flow;
        inFlow.T=outFlow.T;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
                extent={{-100,42},{100,-52}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-100,32},{-100,-42},{100,0},{98,0},{-100,32}},
                lineColor={28,108,200},
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-48,10},{38,-16}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid,
                textString="泵")}), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end Pump_liquid;

      model Tank_liquid

        Interface.SimpleFlange_a inFlow
          annotation (Placement(transformation(extent={{52,32},{72,52}})));
        Interface.SimpleFlange_b outFlow
          annotation (Placement(transformation(extent={{-58,-66},{-38,-46}})));

          //Modelica.SIunits.Length L;
          Modelica.SIunits.Temperature T;
          Modelica.SIunits.MassFlowRate m_flow;
          parameter Modelica.SIunits.Area Area;
          parameter Modelica.SIunits.Density  rou;
          parameter Modelica.SIunits.SpecificHeatCapacity Cp;
          parameter Modelica.SIunits.Length L0=1;
          parameter Modelica.SIunits.Temperature T0;
      initial equation
        //L=L0;
        T=T0;
      equation
        inFlow.m_flow=m_flow;
        outFlow.m_flow=-m_flow;
        //inFlow.m_flow+outFlow.m_flow=rou*Area*der(L);
        Cp*(inFlow.T-outFlow.T)*inFlow.m_flow=rou*Area*L0*Cp*der(T);
       outFlow.T=T;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
                extent={{-52,78},{66,-88}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-52,18},{66,-90}},
                lineColor={28,108,200},
                fillColor={15,55,140},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-28,-50},{46,-74}},
                lineColor={255,255,255},
                fillColor={15,55,140},
                fillPattern=FillPattern.Solid,
                textString="Tank")}), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end Tank_liquid;

      model mixer

        Interface.SimpleFlange_a inFlow1
          annotation (Placement(transformation(extent={{-62,24},{-42,44}})));
        Interface.SimpleFlange_a inFlow2
          annotation (Placement(transformation(extent={{48,24},{68,44}})));
        Interface.SimpleFlange_b outFlow
          annotation (Placement(transformation(extent={{-10,-36},{10,-16}})));

      equation
        inFlow1.m_flow+inFlow2.m_flow=-outFlow.m_flow;
        inFlow1.m_flow*inFlow1.T+inFlow2.m_flow*inFlow2.T=-outFlow.m_flow*outFlow.T;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-98,50},{98,-44}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid), Text(
                extent={{-22,16},{34,-6}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid,
                textString="Mixer")}), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end mixer;

      model sens_liquid

        Interface.SimpleFlange_a inFlow
          annotation (Placement(transformation(extent={{-66,-8},{-46,12}})));
        Interface.SimpleFlange_b outFlow
          annotation (Placement(transformation(extent={{42,-8},{62,12}})));
        Modelica.Blocks.Interfaces.RealOutput y annotation (Placement(transformation(
                extent={{10,38},{30,58}}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-4,32})));

      equation
              outFlow.m_flow+inFlow.m_flow=0;
              outFlow.T=inFlow.T;
              y=inFlow.T;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-80,26},{78,-20}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid), Text(
                extent={{-34,10},{32,-8}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid,
                textString="温度计")}), Diagram(coordinateSystem(
                preserveAspectRatio=false)));
      end sens_liquid;

      model pump_type3
      replaceable package medium =TJUthermo.Media constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model"
        annotation(choicesAllMatching=true);
        Interface.Flange_a InFlow
          annotation (Placement(transformation(extent={{-52,-34},{-32,-14}})));
        Interface.Flange_b OutFlow
          annotation (Placement(transformation(extent={{22,36},{42,56}})));
          parameter Real n_pump(min=0,max=1)=0.67;
          Modelica.SIunits.MassFlowRate m_flow;
          parameter Modelica.SIunits.Pressure p0=0.3e6;
          parameter Modelica.SIunits.Pressure pstart_out=1.2e6;
          parameter Modelica.SIunits.Temperature T0;
          parameter Boolean steadyInFlow=false
          "进口在一段时间内设为稳态";
          parameter Modelica.SIunits.Time steadytime=5
          "开始稳定段设置时间";
          parameter Modelica.SIunits.Volume v_s;
          Modelica.SIunits.Frequency Np;
          Modelica.SIunits.Density rou_in;
          Modelica.SIunits.SpecificEnthalpy hin(start=medium.specificEnthalpy_pT(p0,T0));
          Modelica.SIunits.SpecificEnthalpy hout(start=medium.specificEnthalpy_pT(pstart_out,T0));
          Modelica.SIunits.SpecificEnthalpy hout_s(start=medium.specificEnthalpy_pT(p0,T0));
          //medium.ThermodynamicState state_in;
          //medium.ThermodynamicState state_out;
          Modelica.SIunits.Pressure pin(start=p0);
          Modelica.SIunits.Pressure pout(start=pstart_out);
          Modelica.SIunits.SpecificEntropy s_in;
      //initial equation
        //hin=medium.specificEnthalpy_pT(p0,T0);
        //state_in=medium.setState_pT(p0,T0);
        //state_out=medium.setState_pT(pstart_out,T0);
        Modelica.Blocks.Interfaces.RealInput u annotation (Placement(transformation(
                extent={{32,-8},{72,32}}), iconTransformation(extent={{50,-8},{10,32}})));
      equation
        //InFlow.m_flow=m_flow;
        u=Np;
        m_flow=v_s*rou_in*Np;
        OutFlow.m_flow=-m_flow;
        if steadyInFlow then
          if time<steadytime then
          hin=medium.specificEnthalpy_pT(p0,T0);
          else
        hin=inStream(InFlow.h);
          end if;
        else
        hin=inStream(InFlow.h);
        end if;
        InFlow.h=hin;
        hout=OutFlow.h;
        //OutFlow.h=inStream(OutFlow.h);//修改//
        pin=InFlow.p;
        m_flow=InFlow.m_flow;
        pout=OutFlow.p;
        //state_in=medium.setState_ph(pin,hin);
        rou_in=medium.density_ph(pin,hin);
        //state_out=medium.setState_ph(pout,hout);
        s_in=medium.specificEntropy_ph(pin,hin);
        hout_s=medium.specificEnthalpy_ps(pout,s_in);
        (hout_s-hin)/(hout-hin)=n_pump;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
                extent={{-44,48},{32,-28}},
                lineColor={28,108,200},
                fillColor={85,255,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-8,40},{-44,12},{-8,-20},{-8,-18},{-8,40}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-34,22},{-12,2}},
                lineColor={28,108,200},
                fillColor={85,255,255},
                fillPattern=FillPattern.Solid,
                textString="Pump"),
              Rectangle(
                extent={{-8,30},{50,-6}},
                lineColor={28,108,200},
                fillColor={0,128,255},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-26,-4},{18,-28}},
                lineColor={255,255,255},
                fillColor={0,128,255},
                fillPattern=FillPattern.Solid,
                textString="type 3")}),
                                      Diagram(coordinateSystem(preserveAspectRatio=false)));
      end pump_type3;

      model Turbine_simple_control
      replaceable package medium =TJUthermo.Media constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model"
        annotation(choicesAllMatching=true);
        Interface.Flange_a InFlow
          annotation (Placement(transformation(extent={{-60,32},{-40,52}})));
        Interface.Flange_b OutFlow
          annotation (Placement(transformation(extent={{36,-76},{56,-56}})));
          parameter Real n_turbine(min=0,max=1)=0.67;
          parameter Modelica.SIunits.Volume V_s=0.5;
          Modelica.SIunits.Density rou_in;
          Modelica.SIunits.Frequency Np;
          parameter Modelica.SIunits.Pressure pin_start=1.2e6 annotation(Dialog(tab="Initialization"));
          parameter Modelica.SIunits.Pressure pout_start=0.3e6 annotation(Dialog(tab="Initialization"));
          parameter Modelica.SIunits.Temperature Tin_start=110+273 annotation(Dialog(tab="Initialization"));
         // parameter medium.ThermodynamicState state0=medium.setState_pT(pin_start,Tin_start);
          parameter Modelica.SIunits.SpecificEntropy s_start=medium.specificEntropy_pT(pin_start,Tin_start)
                                                                                                    annotation(Dialog(tab="Initialization"));
          Modelica.SIunits.MassFlowRate m_flow;
          Modelica.SIunits.SpecificEnthalpy hin;//(start=medium.specificEnthalpy_pT(pin_start,Tin_start));
          Modelica.SIunits.SpecificEnthalpy hout;
          Modelica.SIunits.SpecificEnthalpy hout_s;
          medium.ThermodynamicState state_in(phase(
                                             start =     1),s(start=s_start));
                                                            //(start=state0);
          medium.ThermodynamicState state_out;
          Modelica.SIunits.Pressure pin(start=pin_start);
          Modelica.SIunits.Pressure pout;//(start=pout_start);
          Modelica.SIunits.SpecificEntropy s_in(start=s_start);
          Modelica.SIunits.Power W_turbine;
      //initial equation
        //pin=pin_start;
        //hin=medium.specificEnthalpy_pT(pin_start,Tin_start);
          //state_in=medium.setState_pT(pin_start,Tin_start);
          //state_out=medium.setState_ps(pout_start,s_start);
          //hout=medium.specificEnthalpy_ps(pout_start,s_start);
         //hout_s=medium.specificEnthalpy_ps(pout_start,s_start);
        Modelica.Blocks.Interfaces.RealInput u annotation (Placement(transformation(
                extent={{74,-14},{34,26}}), iconTransformation(extent={{56,-6},{30,20}})));
      equation
        u=Np;
        m_flow=V_s*rou_in*Np;
        InFlow.m_flow=m_flow;
        OutFlow.m_flow=-m_flow;
        hin=inStream(InFlow.h);
        hout=OutFlow.h;
        InFlow.h=hin;
        pin=InFlow.p;
        pout=OutFlow.p;
        state_in=medium.setState_ph(pin,hin);
        rou_in=state_in.d;
        state_out=medium.setState_ph(pout,hout);
        s_in=state_in.s;
        hout_s=medium.specificEnthalpy_ps(pout,s_in);
        (hout-hin)/(hout_s-hin)=n_turbine;
        W_turbine=m_flow*(hin-hout);
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Polygon(
                points={{-52,42},{46,88},{46,-68},{-52,-22},{-52,42}},
                lineColor={28,108,200},
                fillColor={85,255,255},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-40,34},{36,26}},
                lineColor={28,108,200},
                fillColor={0,255,0},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-38,30},{32,-14}},
                lineColor={28,108,200},
                fillColor={85,255,255},
                fillPattern=FillPattern.Solid,
                textString="Turbine"),
              Rectangle(
                extent={{-40,-6},{36,-14}},
                lineColor={28,108,200},
                fillColor={0,255,0},
                fillPattern=FillPattern.Solid)}),                      Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={
              Polygon(
                points={{-50,40},{48,86},{48,-70},{-50,-24},{-50,40}},
                lineColor={28,108,200},
                fillColor={85,255,255},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-36,30},{34,-14}},
                lineColor={28,108,200},
                fillColor={85,255,255},
                fillPattern=FillPattern.Solid,
                textString="Turbine"),
              Rectangle(
                extent={{-38,34},{38,26}},
                lineColor={28,108,200},
                fillColor={0,255,0},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-38,-6},{38,-14}},
                lineColor={28,108,200},
                fillColor={0,255,0},
                fillPattern=FillPattern.Solid)}));
      end Turbine_simple_control;

      model Hx1CrossFlow
        replaceable package Medium=ExternalMedia.Media.CoolPropMedium (mediumName="R134A");

         //===============常量================//
           parameter Integer N(min=1)=2 "分割单元体数目";
           parameter Integer np(min=1)=1 "管程数目";

            parameter Modelica.SIunits.CoefficientOfHeatTransfer U=800
          "流体换热系数";
            parameter Modelica.SIunits.Area F_liquid=6.3 "流体换热面积";
            parameter Modelica.SIunits.Volume V_liquid=0.7263 "流体容积";
            parameter Modelica.SIunits.Length L=1;
            parameter Modelica.SIunits.Length dh=0.03;
            parameter Modelica.SIunits.Mass M "流体质量";
            parameter Modelica.SIunits.SpecificHeatCapacity Cp "流体比热";
            parameter Modelica.SIunits.Mass M_wall "管壁质量";
            parameter Modelica.SIunits.SpecificHeatCapacity Cp_wall
          "管壁比热";
            parameter Real max_drhode=100
          "控制密度随压力焓值变化的偏导数值，推荐取100";
            parameter Real witdth_x=0.05 "减缓突变设定值";
            //===============初值区=================//
            parameter Modelica.SIunits.Temperature Tin_start=303 annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.Temperature Tout_start annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.Temperature T1_wall=0.5*(Tout_start+T1_liquid)
          "墙壁初始温度"                                                                            annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.Temperature Tn_wall=0.5*(Tin_start+Tn_liquid)
          "墙壁初始温度"                                                                             annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.Temperature T1_liquid
          "流体入口温度"                                                    annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.Temperature Tn_liquid
          "流体出口温度"                                                    annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.MassFlowRate m_working=0.5
          "工质质量流量"                                                         annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.MassFlowRate m_liquid=0.5
          "流体质量流量"                                                        annotation(Dialog(tab="初始化参数"));
            parameter Modelica.SIunits.Pressure pstart=1.2e6 annotation(Dialog(tab="初始化参数"));
            //Modelica.SIunits.HeatFlowRate Q;
          //===============变量================//
        //接口//
        source_sink_ports.DoubleWalls doubleWalls(N=N,
          final T1=T1_wall,
          final Tn=Tn_wall,
          final M=M_wall,
          final Cp=Cp_wall)
          annotation (Placement(transformation(extent={{-36,-6},{32,22}})));
        Interface.SimpleFlangeMulti LiquidIn(N=N)
          annotation (Placement(transformation(extent={{-10,-92},{10,-72}}),
              iconTransformation(extent={{-10,-92},{10,-72}})));
        Interface.SimpleFlangeMulti LiquidOut(N=N)
          annotation (Placement(transformation(extent={{-10,66},{10,86}}),
              iconTransformation(extent={{-10,66},{10,86}})));
        Interface.Flange_a InFlow
          annotation (Placement(transformation(extent={{-112,-14},{-92,6}}),
              iconTransformation(extent={{-112,-14},{-92,6}})));
        Interface.Flange_b OutFlow
          annotation (Placement(transformation(extent={{92,-12},{112,8}}),
              iconTransformation(extent={{92,-12},{112,8}})));

        units.CrossFlow_liquid crossFlow_liquid(
          final N=N,
          final F=F_liquid/np,
          final V=V_liquid/np,
          final T1=T1_liquid,
          final Tn=Tn_liquid,
          final U=U,
          final Cp=Cp,
          final M=M/np)
          annotation (Placement(transformation(extent={{-52,68},{38,36}})));
        units.veryFlow veryFlow(
          redeclare final package Medium=Medium,
          final N=N,
          final np=np,
          final L=L,
          final dh=dh,
          final max_drdp=max_drhode,
          final max_drdh=max_drhode,
          final witdth_x=witdth_x,
          final Tin_start=Tin_start,
          final Tout_start=Tout_start,
          final m_flow_start=m_working/np,
          final pstart=pstart)
          annotation (Placement(transformation(extent={{-36,-58},{34,-24}})));
      equation
        connect(crossFlow_liquid.node_wall, doubleWalls.Wall_side1) annotation (Line(
              points={{-8.35,47.04},{-8.35,28.52},{-2,28.52},{-2,11.08}}, color={28,108,
                200}));
        connect(LiquidIn, crossFlow_liquid.node_a) annotation (Line(points={{0,-82},{-76,
                -82},{-76,52},{-52,52}},     color={28,108,200}));
        connect(crossFlow_liquid.node_b, LiquidOut) annotation (Line(points={{37.1,51.68},
                {65.55,51.68},{65.55,76},{0,76}},  color={28,108,200}));
        connect(veryFlow.Wall, doubleWalls.Wall_side2) annotation (Line(points={{-0.3,
                -32.84},{-0.3,-12.42},{-2,-12.42},{-2,7.3}}, color={28,108,200}));
        connect(InFlow, veryFlow.InFlow) annotation (Line(points={{-102,-4},{-68,-4},{
                -68,-40.32},{-34.6,-40.32}},  color={0,0,255}));
        connect(veryFlow.OutFlow, OutFlow) annotation (Line(points={{33.3,-39.98},{67.65,
                -39.98},{67.65,-2},{102,-2}},   color={0,0,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
                extent={{-74,36},{80,-46}},
                lineColor={28,108,200},
                fillColor={213,170,255},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-58,48},{-52,-60}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-40,48},{-34,-60}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-20,48},{-14,-60}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{0,48},{6,-60}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{20,48},{26,-60}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{40,48},{46,-60}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{60,48},{66,-60}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-100,2},{-72,-8}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{70,2},{98,-8}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-82,8},{-80,4},{-70,-2},{-80,-14},{-60,-2},{-82,8}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{84,8},{102,-2},{86,-16},{84,8}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-12,-76},{0,-58},{16,-76},{-12,-76}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-12,50},{0,68},{16,50},{-12,50}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-44,24},{54,-32}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid,
                textString="CrossFliud")}),                            Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end Hx1CrossFlow;

      model Pump_type4
      replaceable package medium =TJUthermo.Media constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model"
        annotation(choicesAllMatching=true);
        Interface.Flange_a InFlow
          annotation (Placement(transformation(extent={{-52,-34},{-32,-14}})));
        Interface.Flange_b OutFlow
          annotation (Placement(transformation(extent={{22,36},{42,56}})));
          //parameter Real n_pump(min=0,max=1)=0.67;
          Real eta_pump;
          Real Xpp;
          Real dxp;
          Modelica.SIunits.MassFlowRate m_flow;
          parameter Modelica.SIunits.Pressure p0=0.3e6;
          parameter Modelica.SIunits.Pressure pstart_out=1.2e6;
          parameter Modelica.SIunits.Temperature T0;
          parameter Boolean steadyInFlow=false
          "进口在一段时间内设为稳态";
          parameter Modelica.SIunits.Time steadytime=5
          "开始稳定段设置时间";
          parameter Modelica.SIunits.Volume v_s;
          parameter Modelica.SIunits.Volume v_max;
          Modelica.SIunits.Frequency Np;
          Modelica.SIunits.Density rou_in;
          Modelica.SIunits.SpecificEnthalpy hin(start=medium.specificEnthalpy_pT(p0,T0));
          Modelica.SIunits.SpecificEnthalpy hout(start=medium.specificEnthalpy_pT(pstart_out,T0));
          Modelica.SIunits.SpecificEnthalpy hout_s(start=medium.specificEnthalpy_pT(p0,T0));
          //medium.ThermodynamicState state_in;
          //medium.ThermodynamicState state_out;
          Modelica.SIunits.Pressure pin(start=p0);
          Modelica.SIunits.Pressure pout(start=pstart_out);
          Modelica.SIunits.SpecificEntropy s_in;
          Modelica.SIunits.Power W_pump;
          Modelica.SIunits.Power W_pump_h;
      //initial equation
        //hin=medium.specificEnthalpy_pT(p0,T0);
        //state_in=medium.setState_pT(p0,T0);
        //state_out=medium.setState_pT(pstart_out,T0);
        Modelica.Blocks.Interfaces.RealInput u annotation (Placement(transformation(
                extent={{28,-8},{68,32}}), iconTransformation(extent={{54,-8},{14,32}})));
      equation
        Np=u;
        //InFlow.m_flow=m_flow;
        m_flow=v_s*rou_in*Np;
        OutFlow.m_flow=-m_flow;
        if steadyInFlow then
          if time<steadytime then
          hin=medium.specificEnthalpy_pT(p0,T0);
          else
        hin=inStream(InFlow.h);
          end if;
        else
        hin=inStream(InFlow.h);
        end if;
        InFlow.h=hin;
        hout=OutFlow.h;
        //OutFlow.h=inStream(OutFlow.h);//修改//
        pin=InFlow.p;
        m_flow=InFlow.m_flow;
        pout=OutFlow.p;
        //state_in=medium.setState_ph(pin,hin);
        rou_in=medium.density_ph(pin,hin);
        Xpp=m_flow/(rou_in*v_max);
        if Xpp<0.1 then
          dxp=log(0.1)/log(10);
        elseif Xpp<1 then
          dxp=log(Xpp)/log(10);
        else
          dxp=log(1)/log(10);
        end if;
        eta_pump=0.93-0.11*dxp-0.2*dxp^2-0.06*dxp^3;
        //state_out=medium.setState_ph(pout,hout);
        s_in=medium.specificEntropy_ph(pin,hin);
        hout_s=medium.specificEnthalpy_ps(pout,s_in);
        (hout_s-hin)/(hout-hin)=eta_pump;
        W_pump=Np*v_s*(pout-pin)*eta_pump;
        W_pump_h=m_flow*(hout-hin);
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
                extent={{-44,48},{32,-28}},
                lineColor={28,108,200},
                fillColor={85,255,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-8,40},{-44,12},{-8,-20},{-8,-18},{-8,40}},
                lineColor={28,108,200},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-34,22},{-12,2}},
                lineColor={28,108,200},
                fillColor={85,255,255},
                fillPattern=FillPattern.Solid,
                textString="new
Pump"),       Rectangle(
                extent={{-8,30},{50,-6}},
                lineColor={28,108,200},
                fillColor={0,128,255},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-2,24},{42,0}},
                lineColor={255,255,255},
                fillColor={0,128,255},
                fillPattern=FillPattern.Solid,
                textString="type 3")}),
                                      Diagram(coordinateSystem(preserveAspectRatio=false)));
      end Pump_type4;

      model SolarCollector_New

        Interface.SimpleFlange_a inFlow
          annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
        Interface.SimpleFlange_b outFlow
          annotation (Placement(transformation(extent={{90,-8},{110,12}})));
        Modelica.Blocks.Interfaces.RealInput windspeed annotation (Placement(transformation(
              extent={{-20,-20},{20,20}},
              rotation=-90,
              origin={-78,32}), iconTransformation(
              extent={{-8,-8},{8,8}},
              rotation=-90,
              origin={-72,58})));
        Modelica.Blocks.Interfaces.RealInput I_in annotation (Placement(transformation(
              extent={{-20,-20},{20,20}},
              rotation=-90,
              origin={-18,34}), iconTransformation(
              extent={{-8,-8},{8,8}},
              rotation=-90,
              origin={-12,58})));
        Modelica.Blocks.Interfaces.RealInput Tam_in annotation (Placement(transformation(
              extent={{-20,-20},{20,20}},
              rotation=-90,
              origin={34,70}), iconTransformation(
              extent={{-8,-8},{8,8}},
              rotation=-90,
              origin={46,58})));
          TJUthermo.DataRecord.Arrays Summary(n=N,T=Cell[:].T);
          Modelica.SIunits.HeatFlowRate Q;
          parameter Integer N(min=2) "划分单元体网格数";
          parameter Modelica.SIunits.Area Area "集热面积";
          parameter Real K "吸收系数";
          parameter Modelica.SIunits.Mass M "槽子中导热油质量";
          parameter Real K_cp=0 "导热油比热斜率";
          parameter Real b_cp=2.2e3 "导热油比热截距";
          parameter Real K_rou=0 "导热油比热斜率";
          parameter Real b_rou=900 "导热油比热截距";
          //parameter Modelica.SIunits.SpecificHeatCapacity Cp "导热油比热";
          parameter Modelica.SIunits.MassFlowRate m_flow=1
          "导热油质量流量";
          parameter Modelica.SIunits.Temperature T0 "初始温度";

           parameter Boolean use_eff2=true
          "是否计算所在地入射角偏差";
          parameter Modelica.SIunits.Conversions.NonSIunits.Time_hour ti=10
          "开始模拟太阳时刻";
          parameter Modelica.SIunits.Conversions.NonSIunits.Angle_deg lati=30
          "所在地纬度";
          parameter Integer Date=100 "日期";
        TJUthermo.Components.units.Cell_collector_new Cell[N](
        each A=Area/N,
        each K=K,
        each K_cp=K_cp,
        each b_cp=b_cp,
        each K_rou=K_rou,
        each b_rou=b_rou,
        each M=M/N,
        each T0=T0,
        each use_eff2=use_eff2,
        each ti=ti,
        each lati=lati,
        each Date=Date);

      equation
        if cardinality(windspeed) == 0 then
          windspeed=2.5;
        end if;
         connect(inFlow,Cell[1].inFlow);
         for i in 1:N-1 loop
           connect(Cell[i].outFlow,Cell[i+1].inFlow);
         end for;
         connect(Cell[N].outFlow,outFlow);
         for i in 1:N loop
         connect(I_in,Cell[i].I_in);
         connect(Tam_in,Cell[i].Tam_in);
         //connect(windspeed,Cell[i].windspeed);
         end for;
      algorithm
        Q:=0;
        for i in 1:N loop
        Q:=Q + Cell[i].Q;
        end for;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Polygon(
                points={{-96,58},{-86,48},{-78,32},{-76,20},{-72,-6},{-76,-36},{-76,-38},
                    {-86,-50},{-94,-56},{76,-58},{86,-48},{92,-30},{96,-16},{96,0},{96,
                    16},{90,38},{82,52},{78,58},{-96,58}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-100,18},{100,-16}},
                lineColor={28,108,200},
                fillColor={255,255,0},
                fillPattern=FillPattern.HorizontalCylinder),
              Text(
                extent={{-62,46},{46,24}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={255,255,0},
                textString="New"),
              Text(
                extent={{-56,-20},{52,-48}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={255,255,0},
                textString="PTC"),
              Text(
                extent={{-76,48},{-54,40}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={255,255,0},
                textString="风速"),
              Text(
                extent={{-24,80},{6,68}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={255,255,0},
                textString="辐射"),
              Text(
                extent={{54,54},{78,40}},
                lineColor={28,108,200},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={255,255,0},
                textString="温度")}),    Diagram(coordinateSystem(
                preserveAspectRatio=false), graphics={Rectangle(
                extent={{-100,26},{100,-24}},
                lineColor={28,108,200},
                fillColor={255,255,0},
                fillPattern=FillPattern.Solid)}));
      end SolarCollector_New;

      model Tank_rou
      replaceable package medium =TJUthermo.Media constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model"
        annotation(choicesAllMatching=true);
        Interface.Flange_a InFlow
          annotation (Placement(transformation(extent={{50,30},{70,50}})));
        Interface.Flange_b OutFlow
          annotation (Placement(transformation(extent={{-64,-60},{-44,-40}})));
          Modelica.SIunits.MassFlowRate min;
          Modelica.SIunits.MassFlowRate mout;
          parameter Modelica.SIunits.Pressure p0=0.3e6;
          parameter Modelica.SIunits.Temperature T0=40+273.15;
          parameter Modelica.SIunits.Area Aera=0.5;
          parameter Modelica.SIunits.Length L0=0.5;
          Modelica.SIunits.Length L;
          Modelica.SIunits.Density rou(start=medium.density_pT(p0,T0));
          Modelica.SIunits.SpecificEnthalpy h;
          Modelica.SIunits.SpecificEnthalpy hin;
          medium.ThermodynamicState State_tank;
      initial equation
        L=L0;
        h=medium.specificEnthalpy_pT(p0,T0);
      equation
        min=InFlow.m_flow;
        mout=OutFlow.m_flow;
        InFlow.p=p0;
        OutFlow.p=p0;
        hin=inStream(InFlow.h);
        hin=InFlow.h;
        OutFlow.h=h;
          min+mout=Aera*der(rou*L);
          min*hin+mout*h=Aera*rou*der(h*L);
          rou=medium.density_ph(p0,h);//恒压
          State_tank=medium.setState_ph(p0,h);

        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-52,66},{60,-72}},
                lineColor={28,108,200},
                fillColor={170,213,255},
                fillPattern=FillPattern.VerticalCylinder), Text(
                extent={{-36,30},{40,-36}},
                lineColor={28,108,200},
                fillPattern=FillPattern.VerticalCylinder,
                fillColor={170,170,255},
                textString="Tank"),
              Text(
                extent={{-24,-26},{28,-50}},
                lineColor={28,108,200},
                textString="rou")}),                                   Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end Tank_rou;
    end main_components;

    package Controls

      model PIDcontroller

        Modelica.Blocks.Interfaces.RealInput u
          annotation (Placement(transformation(extent={{-86,-42},{-46,-2}})));
        Modelica.Blocks.Interfaces.RealInput u1
          annotation (Placement(transformation(extent={{-86,10},{-46,50}})));
        Modelica.Blocks.Interfaces.RealOutput y
          annotation (Placement(transformation(extent={{10,-4},{30,16}})));
          Real u2;
          parameter Real u0;
      algorithm
        if u1<50+273 then
          u2:=50+273;
        elseif u1<250+273 then
          u2:=u1;
        else
          u2:=250+273;
        end if;
        if u2>u then
        y:=(u2-u)/50;
        else
          y:=u0;
        end if;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
                extent={{-68,64},{10,-66}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-56,28},{-4,-18}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid,
                textString="controller1"),
              Rectangle(
                extent={{-86,-2},{-46,-44}},
                lineColor={28,108,200},
                fillColor={255,255,85},
                fillPattern=FillPattern.Solid)}), Diagram(coordinateSystem(
                preserveAspectRatio=false)));
      end PIDcontroller;

      model TimeModel

        Modelica.Blocks.Interfaces.RealInput u annotation (Placement(transformation(
                extent={{-100,-18},{-60,22}}), iconTransformation(extent={{-100,-18},{
                  -60,22}})));
        Modelica.Blocks.Interfaces.RealOutput y annotation (Placement(transformation(
                extent={{80,-10},{100,10}}), iconTransformation(extent={{80,-10},{100,
                  10}})));
                  parameter Modelica.SIunits.Time dt;
                  parameter Real umin;
                  parameter Real umax;
      equation
        if time<dt then
          y=0.01;
        else
          y=u;
        end if;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-100,40},{100,-36}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid), Text(
                extent={{-52,18},{50,-12}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid,
                textString="Time")}), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end TimeModel;

      model SuperHeat
      replaceable package medium =TJUthermo.Media constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model"
        annotation(choicesAllMatching=true);
        Interface.Flange_a inFlow
          annotation (Placement(transformation(extent={{-110,-8},{-90,12}})));
        Interface.Flange_b outFlow
          annotation (Placement(transformation(extent={{90,-8},{110,12}})));
       medium.SaturationProperties Sat1;
       Modelica.SIunits.Temperature Tsat;
        Modelica.SIunits.Temperature T;
      Modelica.SIunits.SpecificEnthalpy h;
        Modelica.Blocks.Interfaces.RealOutput y
          annotation (Placement(transformation(extent={{-8,22},{12,42}}),
              iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={2,32})));
      equation
        outFlow.m_flow+inFlow.m_flow=0;
        outFlow.p=inFlow.p;
        h=inStream(inFlow.h);
        inFlow.h=h;
        outFlow.h=h;
        Sat1=medium.setSat_p(inFlow.p);
        Tsat=Sat1.Tsat;
        T=medium.temperature_ph(inFlow.p,h);
        y=T-Tsat;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-98,32},{98,-26}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid), Text(
                extent={{-66,16},{56,-14}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid,
                textString="superHeat")}), Diagram(coordinateSystem(
                preserveAspectRatio=false)));
      end SuperHeat;

      model controlsource

        Modelica.Blocks.Interfaces.RealOutput y
          annotation (Placement(transformation(extent={{32,-10},{52,10}})));
      equation
          if time<4000 then
            y=28;
          //elseif time<21000 then
            //y=45;
          else
            y=35;
          end if;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-50,28},{58,-28}},
                lineColor={28,108,200},
                fillColor={170,213,255},
                fillPattern=FillPattern.Solid)}), Diagram(coordinateSystem(
                preserveAspectRatio=false)));
      end controlsource;

      model sens_PT
      replaceable package medium =TJUthermo.Media constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model"
        annotation(choicesAllMatching=true);
        Interface.Flange_a inFlow
          annotation (Placement(transformation(extent={{-110,-8},{-90,12}})));
        Interface.Flange_b outFlow
          annotation (Placement(transformation(extent={{90,-8},{110,12}})));
        Modelica.SIunits.Temperature T;
      Modelica.SIunits.SpecificEnthalpy h;

        Modelica.Blocks.Interfaces.RealOutput y_T annotation (Placement(
              transformation(extent={{30,48},{50,68}}), iconTransformation(extent={{30,
                  48},{50,68}})));
        Modelica.Blocks.Interfaces.RealOutput y_P annotation (Placement(
              transformation(extent={{-48,44},{-28,64}}), iconTransformation(extent={{
                  -38,46},{-58,66}})));
      equation
        y_T=T;
        y_P=1e-5*outFlow.p;
        outFlow.m_flow+inFlow.m_flow=0;
        outFlow.p=inFlow.p;
        h=inStream(inFlow.h);
        inFlow.h=h;
        outFlow.h=h;
        T=medium.temperature_ph(inFlow.p,h);
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-98,32},{98,-26}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid), Text(
                extent={{-66,16},{56,-14}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid,
                textString="sens_PT"),
              Rectangle(
                extent={{-38,66},{-28,22}},
                lineColor={28,108,200},
                fillColor={127,0,0},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{20,68},{30,24}},
                lineColor={28,108,200},
                fillColor={127,0,0},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-60,44},{-40,8}},
                lineColor={28,108,200},
                fillColor={127,0,0},
                fillPattern=FillPattern.Solid,
                textString="P"),
              Text(
                extent={{32,46},{52,10}},
                lineColor={28,108,200},
                fillColor={127,0,0},
                fillPattern=FillPattern.Solid,
                textString="T")}),         Diagram(coordinateSystem(
                preserveAspectRatio=false)));
      end sens_PT;

      model control_unit1

        Modelica.Blocks.Interfaces.RealInput u
          annotation (Placement(transformation(extent={{-102,-20},{-62,20}})));
        Modelica.Blocks.Interfaces.RealOutput y
          annotation (Placement(transformation(extent={{56,-10},{76,10}})));
          parameter Modelica.SIunits.Time t;
          parameter Real P0;
          parameter Real Y0;
          Real e;
          Real dedt;
          Real sume;
      equation
        e=u-P0;
        if time>t then
          y=Y0+2*e+1.5e3*dedt;//1.5e3
          dedt=der(e);
          //e=der(sume);
        else
          y=Y0;
          dedt=0;
          //sume=0;
        end if;
        der(sume) = e/2500;
        //y = sume+u;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-88,40},{68,-52}},
                lineColor={28,108,200},
                fillColor={127,0,0},
                fillPattern=FillPattern.Solid), Text(
                extent={{-60,16},{28,-28}},
                lineColor={28,108,200},
                fillColor={127,0,0},
                fillPattern=FillPattern.Solid,
                textString="Control")}),                               Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end control_unit1;

      model control_unit2

        Modelica.Blocks.Interfaces.RealInput u
          annotation (Placement(transformation(extent={{-102,-20},{-62,20}})));
        Modelica.Blocks.Interfaces.RealOutput y
          annotation (Placement(transformation(extent={{56,-10},{76,10}})));
          parameter Modelica.SIunits.Time t;
          //parameter Real T0;
         // parameter Real Y0;
          Real e;
          Real dedt;
          //Real sume;
      equation
        e=u;
        if time>t then

          y=1000*dedt;//+1.5e3*dedt;
          dedt=der(e);
          //e=der(sume);
        else
          y=0;
          dedt=0;
          //sume=0;
        end if;
        //der(sume) = e/2500;
        //y = sume+u;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Rectangle(
                extent={{-88,40},{68,-52}},
                lineColor={28,108,200},
                fillColor={127,0,0},
                fillPattern=FillPattern.Solid), Text(
                extent={{-60,16},{28,-28}},
                lineColor={28,108,200},
                fillColor={127,0,0},
                fillPattern=FillPattern.Solid,
                textString="Control2")}),                              Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end control_unit2;
    end Controls;
  end Components;

  package HeatTransferModel
  function constHeatTransfer
    //input
     constant Real pi=3.1415;
     input Modelica.SIunits.CoefficientOfHeatTransfer Ul "液态换热系数";
     input Modelica.SIunits.CoefficientOfHeatTransfer Utp "两相换热系数";
     input Modelica.SIunits.CoefficientOfHeatTransfer Uv;
     input Real witdth_x=0.05 "减缓突变设定值";
     input Real x;
     output Modelica.SIunits.CoefficientOfHeatTransfer U;
  algorithm
  if (x<-witdth_x) then
      U:=Ul;
    elseif (x<witdth_x) then
      U:=Ul+0.5*((Utp-Ul)+(Utp-Ul)*sin(x*2*pi/(4*witdth_x)));
    elseif (x<1-witdth_x) then
      U:=Utp;
    elseif (x<1+witdth_x) then
      U:=Utp-0.5*((Utp-Uv)+(Utp-Uv)*sin((x-1)*2*pi/(4*witdth_x)));
    else
      U:=Uv;
      end if;
  end constHeatTransfer;

    function GnielinskiFcn
      constant Real pi=3.1415;
       input Modelica.SIunits.MassFlowRate mor;
       input Modelica.SIunits.Length dh;
       input Modelica.SIunits.DynamicViscosity V;
       input Real Prf;
       input Modelica.SIunits.ThermalConductivity lamda;
       output Modelica.SIunits.CoefficientOfHeatTransfer U;
    protected
                Real Re;
                Real f;
                Real Nuf;
                Modelica.SIunits.Area A;
    algorithm
      A:=dh*dh*pi/4;
      Re:=mor*dh/(A*V);
      f:=(0.79*log(Re)-1.64)^(-2);
      Nuf:=(f/8)*(Re-1000)*Prf/(1+12.7*sqrt(f/8)*(Prf^(2/3)-1));
      U:=Nuf*lamda/dh;
    end GnielinskiFcn;

    function ShahFcn
       input Real x(min=0,max=1);
       input Modelica.SIunits.Pressure P;
       input Modelica.SIunits.Pressure Pc;
       input Modelica.SIunits.CoefficientOfHeatTransfer Ul;
       output Modelica.SIunits.CoefficientOfHeatTransfer U;
    protected
                Real Pre;
    algorithm
      Pre:=P/Pc;
      U:=Ul*((1-x)^0.8+3.8*x^0.76*(1-x)^0.04/(Pre)^0.38);
    end ShahFcn;

    function GnielinskiFcn_medium
      replaceable package medium=ExternalMedia.Media.CoolPropMedium (mediumName="R134A");
      constant Real pi=3.1415;
       input Modelica.SIunits.Temperature T=300;
       input Modelica.SIunits.Pressure P=300e3;
       input Modelica.SIunits.MassFlowRate mor;
       input Modelica.SIunits.Length dh;
       output Modelica.SIunits.CoefficientOfHeatTransfer U;
    protected
                Real Re;
                Real f;
                Modelica.SIunits.Area A;
                Real Nuf;
                Modelica.SIunits.DynamicViscosity V;
                Real Prf;
                Modelica.SIunits.ThermalConductivity lamda;
                medium.ThermodynamicState state;
    algorithm
      A:=dh*dh*pi/4;
      state:=medium.setState_pT(P,T);
      Prf:=medium.prandtlNumber(state);
      V:=medium.dynamicViscosity(state);
      lamda:=medium.thermalConductivity(state);
      Re:=mor*dh/(A*V);
      f:=(0.79*log(Re)-1.64)^(-2);
      Nuf:=(f/8)*(Re-1000)*Prf/(1+12.7*sqrt(f/8)*(Prf^(2/3)-1));
      U:=Nuf*lamda/dh;
    end GnielinskiFcn_medium;

    function ShahFcn_medium
      replaceable package medium=ExternalMedia.Media.CoolPropMedium (mediumName="R134A");
    input Real x(min=0,max=1);
       input Modelica.SIunits.Pressure P;
       input Modelica.SIunits.CoefficientOfHeatTransfer Ul;
       output Modelica.SIunits.CoefficientOfHeatTransfer U;
    protected
                Real Pre;
      Modelica.SIunits.Pressure Pc;
    algorithm
      Pc:=medium.getCriticalPressure();
      Pre:=P/Pc;
      U:=Ul*((1-x)^0.8+3.8*x^0.76*(1-x)^0.04/(Pre)^0.38);
    end ShahFcn_medium;

    function solar_PTC
      input Modelica.SIunits.Temperature dt;
      input Modelica.SIunits.Time ti;
      input Real t0;
      input Modelica.SIunits.Conversions.NonSIunits.Angle_deg lati;
      input Integer N(min=1,max=365);
      input Modelica.SIunits.HeatFlux I;
      output Real eta_new;
    protected
      constant Real pi=3.1415926;
      Real an1;
      Real delta;//赤纬角
      Modelica.SIunits.Angle time_angle;//时间角
      Modelica.SIunits.Angle latitude;
      Real sin_high;
      Real cos_in;
      Modelica.SIunits.Angle in_angle;//入射角
      Real k1;
    algorithm
        an1:=(360*(284+N)/365)/180*pi;
        delta:=(23.45*sin(an1))/180*pi;
        time_angle:=((ti/3600-12+t0)*15)/180*pi;
        latitude:=lati/180*pi;
        sin_high:=sin(latitude)*sin(delta) + cos(latitude)*cos(delta)*cos(
        time_angle);
        cos_in:=sqrt(sin_high*sin_high + cos(delta)*cos(delta)*sin(time_angle)*sin(
        time_angle));
        in_angle:=acos(cos_in)/pi*180;
        k1:=cos_in - 7e-4*in_angle - 36e-6*in_angle*in_angle;
        eta_new:=0.736*k1*0.95-(29.12*dt+9.52e-7*dt^4)/(409.908*I);
    end solar_PTC;
  end HeatTransferModel;

  package DataRecord

    record Arrays
      parameter Integer n;
      Modelica.SIunits.Temperature[n] T;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Arrays;

    record SummaryClass
       replaceable TJUthermo.DataRecord.Arrays T_profiles(n=N);
       Modelica.SIunits.SpecificEnthalpy[N+1] h_node;
       Modelica.SIunits.SpecificEnthalpy[N] h_cell;
       Modelica.SIunits.MassFlowRate[N+1] M_node;
       Real[N] x;
       Modelica.SIunits.Pressure p;
       parameter Integer N;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end SummaryClass;

    type GasConstant = Real (final quantity="GasConstant", unit="J/(mol.K)");
    record SatProperty
      Modelica.SIunits.Pressure Psat;
      Modelica.SIunits.Density roul;
      Modelica.SIunits.Density rouv;
      Modelica.SIunits.SpecificEnthalpy hl;
      Modelica.SIunits.SpecificEnthalpy hv;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end SatProperty;
  end DataRecord;

  package Example
    package Test
      model Test_Cell
        Components.source_sink_ports.InFlow inFlow(
          h0=2.5e5,
          UseT=true,
          redeclare package medium = ThermoCycle.Media.R410a,
          m0=-1.5,
          p0=800000,
          T=283.15)
          annotation (Placement(transformation(extent={{-80,-6},{-60,14}})));
        Components.source_sink_ports.OutFlow outFlow(m0=0.5,
          redeclare package medium = ThermoCycle.Media.R410a,
          p0=1000000)
          annotation (Placement(transformation(extent={{72,-4},{52,16}})));
        Components.source_sink_ports.Wall wall(T=363.15)
          annotation (Placement(transformation(extent={{-28,20},{16,64}})));
        Components.units.Cell cell1(
          F=0.5,
          V=0.06,
          max_drdp=100,
          max_drdh=100,
          redeclare package medium = ThermoCycle.Media.R410a,
          Tin_start=303.15,
          Tout_start=353.15,
          pstart=1000000)
          annotation (Placement(transformation(extent={{-40,-16},{28,22}})));
      equation
        connect(inFlow.node, cell1.A) annotation (Line(points={{-77.6,3.6},{
                -55.9,3.6},{-55.9,3},{-39.32,3}}, color={0,0,255}));
        connect(cell1.B, outFlow.node) annotation (Line(points={{27.32,3},{
                39.66,3},{39.66,5.8},{53.8,5.8}}, color={0,0,255}));
        connect(wall.thermoT, cell1.wall) annotation (Line(points={{-6.22,41.78},
                {-6.22,24.89},{-6.34,24.89},{-6.34,6.99}}, color={28,108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end Test_Cell;

      model TestSimpleFlow
        Components.units.SimpleFlow simpleFlow(
          N=10,
          F=0.5*10,
          V=0.06*10,
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          m_flow_start=0.4,
          Ul=1500,
          Uv=1000,
          Tin_start=358.15,
          Tout_start=373.15)
          annotation (Placement(transformation(extent={{-38,-36},{34,28}})));
        Components.source_sink_ports.InFlow inFlow(
          m0=-0.5,
          UseT=true,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0=1200000,
          T=358.15)
          annotation (Placement(transformation(extent={{-88,-10},{-68,10}})));
        Components.source_sink_ports.OutFlow outFlow(m0=0.5,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0=1200000)
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=0,
              origin={68,0})));
        Components.source_sink_ports.Walls walls(
          T1=368.15,
          Tn=403.15,
          N=10)
          annotation (Placement(transformation(extent={{-36,24},{34,52}})));
      equation
        connect(inFlow.node, simpleFlow.InFlow) annotation (Line(points={{-69.8,
                -0.2},{-53.9,-0.2},{-53.9,-2.72},{-36.56,-2.72}}, color={0,0,
                255}));
        connect(simpleFlow.OutFlow, outFlow.node) annotation (Line(points={{
                33.28,-2.08},{46.64,-2.08},{46.64,-0.2},{59.8,-0.2}}, color={0,
                0,255}));
        connect(walls.walls, simpleFlow.Wall) annotation (Line(points={{-2.05,
                38.28},{-1.28,38.28},{-1.28,11.36}}, color={28,108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestSimpleFlow;

      model TestDoubleFlow
        Components.units.SimpleFlow working_Flow(
          N=10,
          F=5,
          V=0.6,
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          max_drdp=100,
          max_drdh=100,
          witdth_x=0.1,
          Utp=1000,
          Uv=1000,
          Tin_start=358.15,
          Tout_start=381.15,
          pstart=1200000)
          annotation (Placement(transformation(extent={{-38,-68},{46,-28}})));
        Components.source_sink_ports.InFlow inFlow(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          m0=-0.5,
          UseT=true,
          p0=1200000,
          T=358.15)
          annotation (Placement(transformation(extent={{-88,-58},{-68,-38}})));
        Components.source_sink_ports.OutFlow outFlow1(
          m0=0.5,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0=1200000)
          annotation (Placement(transformation(extent={{94,-58},{74,-38}})));
        Components.units.SimpleFlow water_Flow(
          N=10,
          Ul=700,
          F=5,
          V=0.7,
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Tin_start=413.15,
          Tout_start=383.15,
          pstart=1700000)
          annotation (Placement(transformation(extent={{-44,50},{44,10}})));
        Components.source_sink_ports.InFlow inFlow1(
          m0=0.5,
          UseT=true,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0=1700000,
          T=413.15)                                 annotation (Placement(
              transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={-78,26})));
        Components.source_sink_ports.OutFlow outFlow2(
          m0=0.5,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0=1700000)
          annotation (Placement(transformation(extent={{86,16},{66,36}})));
        Components.source_sink_ports.DoubleWalls doubleWalls(
          N=10,
          T1=373.15,
          Tn=383.15)
          annotation (Placement(transformation(extent={{-8,-16},{12,4}})));
      equation
        connect(inFlow.node, working_Flow.InFlow) annotation (Line(points={{-69.8,
                -48.2},{-53.9,-48.2},{-53.9,-47.2},{-36.32,-47.2}},       color=
               {0,0,255}));
        connect(outFlow1.node, working_Flow.OutFlow) annotation (Line(points={{75.8,
                -48.2},{52.9,-48.2},{52.9,-46.8},{45.16,-46.8}},      color={0,
                0,255}));
        connect(inFlow1.node,water_Flow. InFlow) annotation (Line(points={{-69.8,
                26.2},{-69.8,27.1},{-42.24,27.1},{-42.24,29.2}},       color={0,
                0,255}));
        connect(outFlow2.node,water_Flow. OutFlow) annotation (Line(points={{67.8,
                25.8},{64.1,25.8},{64.1,28.8},{43.12,28.8}},          color={0,
                0,255}));
        connect(water_Flow.Wall, doubleWalls.Wall_side1) annotation (Line(
              points={{0.88,20.4},{0.88,8.2},{2,8.2},{2,-3.8}}, color={28,108,
                200}));
        connect(doubleWalls.Wall_side2, working_Flow.Wall) annotation (Line(
              points={{2,-6.5},{4,-6.5},{4,-38.4},{4.84,-38.4}}, color={28,108,
                200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-28,-66},{36,-90}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid,
                textString=
                    "双流式换热器（问题：不能区分正逆流）"), Text(
                extent={{-40,74},{46,48}},
                lineColor={28,108,200},
                textString="仿真有问题")}));
      end TestDoubleFlow;

      model Test_workingfluids
      replaceable package medium = ExternalMedia.Media.CoolPropMedium (
        mediumName="R410a",substanceNames={"REFPROP-MIX:R32[0.697615]&R125[0.302385]"}) annotation(choicesFromPackage=true);
        parameter String mediumName="R123";
        parameter Modelica.SIunits.Temperature T=319;
            parameter Modelica.SIunits.Pressure p=0.6e6;
       medium.ThermodynamicState state;
      equation
        state=medium.setState_pT(p,T)
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end Test_workingfluids;

      model TestSimpleFlow_cooling
        Components.units.SimpleFlow simpleFlow(
          N=10,
          F=0.5*10,
          V=0.06*10,
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          m_flow_start=0.4,
          Ul=1500,
          Uv=1000,
          Tin_start=363.15,
          Tout_start=348.15)
          annotation (Placement(transformation(extent={{-28,-26},{44,38}})));
        Components.source_sink_ports.InFlow inFlow(
          m0=-0.5,
          UseT=true,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0=1200000,
          T=358.15)
          annotation (Placement(transformation(extent={{-78,0},{-58,20}})));
        Components.source_sink_ports.OutFlow outFlow(m0=0.5,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0=1200000)
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=0,
              origin={78,10})));
        Components.source_sink_ports.Walls walls(
          N=10,
          T1=358.15,
          Tn=343.15)
          annotation (Placement(transformation(extent={{-26,34},{44,62}})));
      equation
        connect(inFlow.node,simpleFlow. InFlow) annotation (Line(points={{-59.8,
                9.8},{-43.9,9.8},{-43.9,7.28},{-26.56,7.28}},     color={0,0,
                255}));
        connect(simpleFlow.OutFlow,outFlow. node) annotation (Line(points={{43.28,
                7.92},{56.64,7.92},{56.64,9.8},{69.8,9.8}},           color={0,
                0,255}));
        connect(walls.walls,simpleFlow. Wall) annotation (Line(points={{7.95,
                48.28},{8.72,48.28},{8.72,21.36}},   color={28,108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestSimpleFlow_cooling;

      model Test_Cellliquid
        Components.source_sink_ports.Wall wall(T=293.15)
          annotation (Placement(transformation(extent={{-8,32},{12,52}})));
        Components.source_sink_ports.InFlow_Tmflow inFlow_Tmflow(
          m0=-1,
          T0=323.15,
          use_m0=true)
          annotation (Placement(transformation(extent={{-74,-32},{-54,-12}})));
        Components.source_sink_ports.OutFlow_Tmflow outFlow_Tmflow(m0=1, use_m0=
             false)
          annotation (Placement(transformation(extent={{62,-38},{82,-18}})));
        Components.units.Cell_liquid cell_liquid1(
          F=0.5,
          V=0.06,
          U=700,
          Cp=2.1e3,
          M=1)
          annotation (Placement(transformation(extent={{-20,-14},{22,14}})));
      equation
        connect(wall.thermoT, cell_liquid1.wall) annotation (Line(points={{1.9,
                41.9},{1.9,20.95},{0.79,20.95},{0.79,2.94}}, color={28,108,200}));
        connect(inFlow_Tmflow.outFlow, cell_liquid1.A) annotation (Line(points=
                {{-54,-22},{-34,-22},{-34,1.77636e-015},{-19.58,1.77636e-015}},
              color={28,108,200}));
        connect(cell_liquid1.B, outFlow_Tmflow.inFlow) annotation (Line(points=
                {{21.58,1.77636e-015},{44,1.77636e-015},{44,-28.2},{62,-28.2}},
              color={28,108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end Test_Cellliquid;

      model Test_twoCell
        Components.units.Cell cell(
          F=0.5,
          V=0.06,
          hstart=3.3e5,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          Tstart=353.15)
          annotation (Placement(transformation(extent={{-46,-72},{14,-24}})));
        Components.source_sink_ports.InFlow inFlow(
          m0=-0.5,
          h0=2.5e5,
          UseT=true,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0=1200000,
          T=363.15)
          annotation (Placement(transformation(extent={{-90,-58},{-70,-38}})));
        Components.source_sink_ports.OutFlow outFlow(m0=0.5,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0=1200000)
          annotation (Placement(transformation(extent={{130,-60},{110,-40}})));
        Components.units.Cell_liquid cell_liquid(
          F=0.5,
          V=0.06,
          Cp=2.1e3,
          m_flow=0.4,
          M=1,
          Tstart=423.15,
          U=1500)
          annotation (Placement(transformation(extent={{-32,24},{14,54}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid(T0=423.15)
          annotation (Placement(transformation(extent={{-84,26},{-64,46}})));
        Components.units.Cell cell1(
          F=0.5,
          V=0.06,
          hstart=3.3e5,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          Tstart=353.15)
          annotation (Placement(transformation(extent={{34,-72},{94,-24}})));
        Components.units.Cell_liquid cell_liquid1(
          F=0.5,
          V=0.06,
          Cp=2.1e3,
          m_flow=0.4,
          M=1,
          Tstart=423.15,
          U=1500)
          annotation (Placement(transformation(extent={{30,26},{76,56}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid(m0=1)
          annotation (Placement(transformation(extent={{100,30},{120,50}})));
        Components.source_sink_ports.Double_wall double_wall(
          Cp=1.4e3,
          T0=383.15,
          M=0.2)
          annotation (Placement(transformation(extent={{-28,-16},{-8,4}})));
        Components.source_sink_ports.Double_wall double_wall1(
          Cp=1.4e3,
          T0=383.15,
          M=0.2)
          annotation (Placement(transformation(extent={{48,-18},{68,2}})));
      equation
        connect(inFlow.node,cell. A) annotation (Line(points={{-71.8,-48.2},{
                -46.9,-48.2},{-46.9,-48},{-45.4,-48}},
                                                 color={0,0,255}));
        connect(inFlow_liquid.A, cell_liquid.A) annotation (Line(points={{-64.6,
                35.8},{-43.3,35.8},{-43.3,39},{-31.54,39}}, color={28,108,200}));
        connect(cell_liquid.B, cell_liquid1.A) annotation (Line(points={{13.54,
                39},{21.77,39},{21.77,41},{30.46,41}}, color={28,108,200}));
        connect(cell.B, cell1.A) annotation (Line(points={{13.4,-48},{26,-48},{
                34.6,-48}}, color={0,0,255}));
        connect(cell1.B, outFlow.node) annotation (Line(points={{93.4,-48},{
                93.4,-47},{111.8,-47},{111.8,-50.2}}, color={0,0,255}));
        connect(cell_liquid1.B, outFlow_liquid.A) annotation (Line(points={{
                75.54,41},{87.77,41},{87.77,40},{100,40}}, color={28,108,200}));
        connect(cell_liquid.wall, double_wall.A) annotation (Line(points={{
                -9.23,42.15},{-9.23,20.075},{-18,20.075},{-18,-3.6}}, color={28,
                108,200}));
        connect(double_wall.B, cell.wall) annotation (Line(points={{-18,-9.4},{
                -14,-9.4},{-14,-42.96},{-16.3,-42.96}}, color={28,108,200}));
        connect(cell_liquid1.wall, double_wall1.A) annotation (Line(points={{
                52.77,44.15},{52.77,19.075},{58,19.075},{58,-5.6}}, color={28,
                108,200}));
        connect(double_wall1.B, cell1.wall) annotation (Line(points={{58,-11.4},
                {62,-11.4},{62,-42.96},{63.7,-42.96}}, color={28,108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end Test_twoCell;

      model TestTwoliquid
        Components.units.Cell_liquid cell_liquid(
          F=0.5,
          V=0.06,
          Cp=2.1e3,
          m_flow=0.4,
          M=1,
          Tstart=423.15,
          U=1500)
          annotation (Placement(transformation(extent={{-52,2},{-6,32}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid(T0=423.15)
          annotation (Placement(transformation(extent={{-104,4},{-84,24}})));
        Components.units.Cell_liquid cell_liquid1(
          F=0.5,
          V=0.06,
          Cp=2.1e3,
          m_flow=0.4,
          M=1,
          Tstart=423.15,
          U=1500)
          annotation (Placement(transformation(extent={{10,4},{56,34}})));
        Components.source_sink_ports.Wall wall2(T=388.15)
          annotation (Placement(transformation(extent={{-30,-40},{-10,-20}})));
        Components.source_sink_ports.Wall wall3(T=388.15)
          annotation (Placement(transformation(extent={{24,-30},{44,-10}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid(m0=1)
          annotation (Placement(transformation(extent={{76,8},{96,28}})));
      equation
        connect(inFlow_liquid.A, cell_liquid.A) annotation (Line(points={{-84.6,
                13.8},{-63.3,13.8},{-63.3,17},{-51.54,17}}, color={28,108,200}));
        connect(cell_liquid.B, cell_liquid1.A) annotation (Line(points={{-6.46,
                17},{1.77,17},{1.77,19},{10.46,19}}, color={28,108,200}));
        connect(wall2.thermoT, cell_liquid.wall) annotation (Line(points={{
                -20.1,-30.1},{-20.1,20.15},{-29.23,20.15}}, color={28,108,200}));
        connect(wall3.thermoT, cell_liquid1.wall) annotation (Line(points={{
                33.9,-20.1},{33.9,0.95},{32.77,0.95},{32.77,22.15}}, color={28,
                108,200}));
        connect(cell_liquid1.B, outFlow_liquid.A) annotation (Line(points={{
                55.54,19},{66.77,19},{66.77,18},{76,18}}, color={28,108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestTwoliquid;

      model TestDoubleFlow2
        Components.units.SimpleFlow working_Flow(
          N=10,
          F=5,
          V=0.6,
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          max_drdp=100,
          max_drdh=100,
          witdth_x=0.1,
          Utp=1000,
          Uv=1000,
          Tin_start=358.15,
          Tout_start=381.15,
          pstart=1200000)
          annotation (Placement(transformation(extent={{-30,-58},{54,-18}})));
        Components.source_sink_ports.InFlow inFlow(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          m0=-0.5,
          UseT=true,
          p0=1200000,
          T=358.15)
          annotation (Placement(transformation(extent={{-78,-48},{-58,-28}})));
        Components.source_sink_ports.OutFlow outFlow1(
          m0=0.5,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0=1200000)
          annotation (Placement(transformation(extent={{104,-48},{84,-28}})));
        Components.units.SimpleFlowliquid simpleFlowliquid(
          N=10,
          U=800,
          Cp=2.1e3,
          M=0.1,
          F=5,
          V=0.6,
          m_flow=0.5,
          T1=413.15,
          Tn=383.15)
          annotation (Placement(transformation(extent={{-36,-2},{56,54}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid(T0=413.15)
          annotation (Placement(transformation(extent={{-84,0},{-64,20}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid(m0=1)
          annotation (Placement(transformation(extent={{90,6},{110,26}})));
        Components.source_sink_ports.DoubleWalls doubleWalls(
          N=10,
          T1=383.15,
          Tn=383.15,
          M=1,
          Cp=2e3)
          annotation (Placement(transformation(extent={{-2,-14},{18,6}})));
      equation
        connect(inFlow.node,working_Flow. InFlow) annotation (Line(points={{-59.8,
                -38.2},{-43.9,-38.2},{-43.9,-37.2},{-28.32,-37.2}},       color=
               {0,0,255}));
        connect(outFlow1.node,working_Flow. OutFlow) annotation (Line(points={{85.8,
                -38.2},{62.9,-38.2},{62.9,-36.8},{53.16,-36.8}},      color={0,
                0,255}));
        connect(inFlow_liquid.A,simpleFlowliquid. node_a) annotation (Line(
              points={{-64.6,9.8},{-47.3,9.8},{-47.3,26},{-36,26}}, color={28,
                108,200}));
        connect(simpleFlowliquid.node_b,outFlow_liquid. A) annotation (Line(
              points={{55.08,26.56},{76.54,26.56},{76.54,16},{90,16}},
                                                                   color={28,
                108,200}));
        connect(simpleFlowliquid.node_wall, doubleWalls.Wall_side1) annotation (
           Line(points={{8.62,34.68},{8.62,-1.8},{8,-1.8}}, color={28,108,200}));
        connect(doubleWalls.Wall_side2, working_Flow.Wall) annotation (Line(
              points={{8,-4.5},{12,-4.5},{12,-28.4},{12.84,-28.4}}, color={28,
                108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-18,-56},{46,-80}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid,
                textString="双流式换热器（顺流）")}));
      end TestDoubleFlow2;

      model TestFlow_liquids
        Components.units.SimpleFlowliquid simpleFlowliquid(
          N=10,
          F=0.5,
          V=0.06,
          U=800,
          Cp=2.1e3,
          M=0.1,
          T1=353.15,
          Tn=338.15)
          annotation (Placement(transformation(extent={{-42,-26},{50,30}})));
        Components.source_sink_ports.Walls walls(
          N=10,
          T1=328.15,
          Tn=308.15)
          annotation (Placement(transformation(extent={{-10,42},{10,62}})));
        Components.source_sink_ports.InFlow_Tmflow inFlow_Tmflow(
          use_m0=true,
          T0=353.15,
          m0=-0.5)
          annotation (Placement(transformation(extent={{-82,-46},{-62,-26}})));
        Components.source_sink_ports.OutFlow_Tmflow outFlow_Tmflow(m0=1, use_m0=
             false)
          annotation (Placement(transformation(extent={{80,-52},{100,-32}})));
      equation
        connect(walls.walls, simpleFlowliquid.node_wall) annotation (Line(
              points={{-0.3,52.2},{-0.3,32.1},{2.62,32.1},{2.62,10.68}}, color=
                {28,108,200}));
        connect(inFlow_Tmflow.outFlow, simpleFlowliquid.node_a) annotation (
            Line(points={{-62,-36},{-50,-36},{-50,2},{-42,2}}, color={28,108,
                200}));
        connect(simpleFlowliquid.node_b, outFlow_Tmflow.inFlow) annotation (
            Line(points={{49.08,2.56},{49.08,-19.72},{80,-19.72},{80,-42.2}},
              color={28,108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestFlow_liquids;

      model TestCounterFlow
        Components.units.SimpleFlowliquid simpleFlowliquid(
          N=10,
          m_flow=0.5,
          F=6.3,
          V=0.7263,
          U=700,
          Cp=4.2e3,
          M=700,
          T1=423.15,
          Tn=373.15)
          annotation (Placement(transformation(extent={{-52,36},{40,92}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid1(T0=423.15)
          annotation (Placement(transformation(extent={{-94,20},{-74,40}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(m0=0.5)
          annotation (Placement(transformation(extent={{68,12},{88,32}})));
        Components.units.SimpleFlow simpleFlow(
          N=10,
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Uv=700,
          F=6.3,
          V=0.7263,
          max_drdp=200,
          max_drdh=200,
          Tin_start=363.15,
          Tout_start=388.15,
          m_flow_start=0.5)
          annotation (Placement(transformation(extent={{-40,-100},{32,-36}})));
        Components.source_sink_ports.InFlow inFlow2(
          m0=-0.5,
          UseT=true,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0=1200000,
          T=363.15)
          annotation (Placement(transformation(extent={{-84,-92},{-64,-72}})));
        Components.source_sink_ports.OutFlow outFlow(m0=0.5,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0=1200000)
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=0,
              origin={82,-70})));
        Components.source_sink_ports.DoubleWalls doubleWalls(
          N=10,
          T1=388.15,
          Tn=388.15,
          M=20,
          Cp=0.3e3)
          annotation (Placement(transformation(extent={{-18,-8},{2,12}})));
      equation
        connect(inFlow_liquid1.A, simpleFlowliquid.node_a) annotation (Line(
              points={{-74.6,29.8},{-74.6,45.9},{-52,45.9},{-52,64}}, color={28,
                108,200}));
        connect(simpleFlowliquid.node_b, outFlow_liquid1.A) annotation (Line(
              points={{39.08,64.56},{56.54,64.56},{56.54,22},{68,22}}, color={
                28,108,200}));
        connect(inFlow2.node, simpleFlow.InFlow) annotation (Line(points={{
                -65.8,-82.2},{-49.9,-82.2},{-49.9,-66.72},{-38.56,-66.72}},
              color={0,0,255}));
        connect(simpleFlow.OutFlow,outFlow. node) annotation (Line(points={{31.28,
                -66.08},{50.64,-66.08},{50.64,-70.2},{73.8,-70.2}},   color={0,
                0,255}));
        connect(simpleFlowliquid.node_wall, doubleWalls.Wall_side1) annotation (
           Line(points={{-7.38,72.68},{-7.38,38.34},{-8,38.34},{-8,4.2}}, color=
               {28,108,200}));
        connect(doubleWalls.Wall_side2, simpleFlow.Wall) annotation (Line(
              points={{-8,1.5},{-2,1.5},{-2,-52.64},{-3.28,-52.64}}, color={28,
                108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-86,-8},{-20,-36}},
                lineColor={28,108,200},
                textString="顺流式换热器")}));
      end TestCounterFlow;

      model TestCellCollector
        Modelica.Blocks.Sources.Constant const1(k=700)
          annotation (Placement(transformation(extent={{14,60},{34,80}})));
        Modelica.Blocks.Sources.Constant const2(k=25)
          annotation (Placement(transformation(extent={{-24,60},{-4,80}})));
        Components.units.Cell_collector cell_collector1(
          A=15,
          K=1,
          M=1,
          Cp=2.1e3,
          m_flow=0.5,
          T0=343.15)
          annotation (Placement(transformation(extent={{-32,-24},{46,32}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid(T0=343.15)
          annotation (Placement(transformation(extent={{-82,-36},{-62,-16}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid(m0=0.5)
          annotation (Placement(transformation(extent={{78,-44},{98,-24}})));
      equation
        connect(const2.y, cell_collector1.Tam_in) annotation (Line(points={{-3,
                70},{2,70},{2,11.84},{6.22,11.84}}, color={0,0,127}));
        connect(const1.y, cell_collector1.I_in) annotation (Line(points={{35,70},
                {32,70},{32,11},{29.23,11}}, color={0,0,127}));
        connect(inFlow_liquid.A, cell_collector1.inFlow) annotation (Line(
              points={{-62.6,-26.2},{-45.3,-26.2},{-45.3,5.12},{-31.22,5.12}},
              color={28,108,200}));
        connect(cell_collector1.outFlow, outFlow_liquid.A) annotation (Line(
              points={{45.22,5.12},{58.61,5.12},{58.61,-34},{78,-34}}, color={
                28,108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestCellCollector;

      model TestSource
        Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(
          tableOnFile=true,
          tableName="uaua",
          fileName="uaua.mat",
          columns={2})
          annotation (Placement(transformation(extent={{-24,26},{-4,46}})));
        Modelica.Blocks.Sources.CombiTimeTable combiTimeTable1(
          tableOnFile=true,
          columns={2},
          tableName="haha",
          fileName="haha.mat")
          annotation (Placement(transformation(extent={{-26,-32},{-6,-12}})));
        Modelica.Blocks.Sources.CombiTimeTable combiTimeTable2(
          tableOnFile=true,
          columns={2},
          tableName="ua",
          fileName="ua.mat")
          annotation (Placement(transformation(extent={{-26,-72},{-6,-52}})));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestSource;

      model TestGnielinski
        replaceable package medium=ExternalMedia.Media.CoolPropMedium (mediumName="R245fa");
        parameter Modelica.SIunits.Temperature T=300;
        parameter Modelica.SIunits.Pressure P=300e3;
        parameter Modelica.SIunits.MassFlowRate mor=0.8;
        parameter Modelica.SIunits.Length dh=0.03;
        Modelica.SIunits.DynamicViscosity V;
        Modelica.SIunits.Pressure Pc;
        //Modelica.SIunits.DynamicViscosity Vl;
        Real Prf;
        Modelica.SIunits.ThermalConductivity lamda;
        Modelica.SIunits.CoefficientOfHeatTransfer U;
          Modelica.SIunits.CoefficientOfHeatTransfer U1;
            Modelica.SIunits.CoefficientOfHeatTransfer U2;
        medium.ThermodynamicState state;
        medium.SaturationProperties Sat1;
       // medium.ThermodynamicState statel;
      equation
        Prf=medium.prandtlNumber(state);
        V=medium.dynamicViscosity(state);
        lamda=medium.thermalConductivity(state);
        state=medium.setState_pT(P,T);
        Sat1=medium.setSat_p(P);
        Pc=medium.getCriticalPressure();
        //statel=medium.setState_pT(P,Sat1.Tsat-0.5);
        U=HeatTransferModel.GnielinskiFcn(mor=mor,dh=dh,V=V,Prf=Prf,lamda=lamda);
        U1=HeatTransferModel.GnielinskiFcn_medium(mor=mor,dh=dh,T=Sat1.Tsat-0.5,P=P);
        U2=HeatTransferModel.ShahFcn(x=0.5,P=P,Pc=Pc,Ul=U1);
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestGnielinski;

      model Test_varyCell
        Components.units.Cell_varyHTC cell_varyHTC(L=10, dh=0.03)
          annotation (Placement(transformation(extent={{-32,-26},{32,30}})));
        Components.source_sink_ports.InFlow inFlow(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0(displayUnit="kPa") = 500000,
          UseT=true,
          m0=-0.5,
          T=313.15)
          annotation (Placement(transformation(extent={{-88,-4},{-68,16}})));
        Components.source_sink_ports.OutFlow outFlow(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          m0=0.5,
          p0(displayUnit="kPa") = 500000)
          annotation (Placement(transformation(extent={{58,18},{78,38}})));
        Components.source_sink_ports.Wall wall(T=363.15)
          annotation (Placement(transformation(extent={{-36,26},{8,70}})));
      equation
        connect(inFlow.node, cell_varyHTC.A) annotation (Line(points={{-69.8,
                5.8},{-51.9,5.8},{-51.9,2},{-31.36,2}}, color={0,0,255}));
        connect(outFlow.node, cell_varyHTC.B) annotation (Line(points={{76.2,
                27.8},{76.2,16.9},{31.36,16.9},{31.36,2}}, color={0,0,255}));
        connect(wall.thermoT, cell_varyHTC.wall) annotation (Line(points={{
                -14.22,47.78},{-14.22,27.89},{-0.32,27.89},{-0.32,7.88}}, color=
               {28,108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end Test_varyCell;

      model Test_solar_PTC
        parameter Modelica.SIunits.Conversions.NonSIunits.Angle_deg lati=45;
        parameter Modelica.SIunits.Temperature dt=100;
        Real c;
      equation
        c=TJUthermo.HeatTransferModel.solar_PTC(dt=dt,ti=time,t0=12,lati=lati,N=1,I=300)
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));

      end Test_solar_PTC;
    end Test;

    package TestComponets
      model TestHx1counter_liquid
        Components.main_components.Hx1counter_liquid hx1counter_liquid(
          N=10,
          Cp1=4.2e3,
          M1=10,
          M2=10,
          M_wall=1,
          Cp_wall=899,
          Cp2=2.3e3,
          U1=2000,
          U2=1500,
          T1_liquid=303.15,
          T1n_liquid=353.15,
          T2_liquid=393.15,
          T2n_liquid=318.15)
          annotation (Placement(transformation(extent={{-42,-16},{36,40}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid(T0=303.15)
          annotation (Placement(transformation(extent={{-88,-10},{-68,10}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid1(T0=393.15)
          annotation (Placement(transformation(extent={{34,44},{54,64}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid(m0=2)
          annotation (Placement(transformation(extent={{62,-10},{82,10}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(m0=1)
          annotation (Placement(transformation(extent={{-62,50},{-42,70}})));
      equation
        connect(inFlow_liquid.A, hx1counter_liquid.inflow1) annotation (Line(
              points={{-68.6,-0.2},{-54.3,-0.2},{-54.3,0.8},{-42.78,0.8}},
              color={28,108,200}));
        connect(inFlow_liquid1.A, hx1counter_liquid.inflow2) annotation (Line(
              points={{53.4,53.8},{53.4,25.9},{35.22,25.9},{35.22,25.44}},
              color={28,108,200}));
        connect(hx1counter_liquid.outflow1, outFlow_liquid.A) annotation (Line(
              points={{35.22,0.8},{48.61,0.8},{48.61,0},{62,0}}, color={28,108,
                200}));
        connect(outFlow_liquid1.A, hx1counter_liquid.outflow2) annotation (Line(
              points={{-62,60},{-80,60},{-80,26},{-42,26}}, color={28,108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestHx1counter_liquid;

      model TestHx1Simple
        Components.source_sink_ports.InFlow inFlow(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          UseT=true,
          m0=-0.5,
          p0=1200000,
          T=363.15)
          annotation (Placement(transformation(extent={{-72,-36},{-52,-16}})));
        Components.source_sink_ports.OutFlow outFlow(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          m0=0.5,
          p0=1200000)
          annotation (Placement(transformation(extent={{46,-46},{66,-26}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid(T0=423.15)
          annotation (Placement(transformation(extent={{-82,10},{-62,30}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid(m0=0.5)
          annotation (Placement(transformation(extent={{60,14},{80,34}})));
        Components.main_components.Hx1simple hx1simple1(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          N=10,
          Uv=700,
          U=700,
          M=0.7e3,
          Cp=4.2e3,
          M_wall=20,
          Cp_wall=0.3e3,
          max_drhode=200,
          m_liquid=0.5,
          Tin_start=363.15,
          Tout_start=388.15,
          T1_liquid=423.15,
          Tn_liquid=403.15)
          annotation (Placement(transformation(extent={{-28,-12},{20,30}})));
      equation
        connect(inFlow.node, hx1simple1.InFlow_working) annotation (Line(points={{-53.8,
                -26.2},{-36.9,-26.2},{-36.9,2.7},{-28,2.7}},        color={0,0,
                255}));
        connect(hx1simple1.OutFlow_working, outFlow.node) annotation (Line(
              points={{19.52,2.7},{42.76,2.7},{42.76,-36.2},{64.2,-36.2}},
              color={0,0,255}));
        connect(outFlow_liquid.A, hx1simple1.OutFlow_liquid) annotation (Line(
              points={{60,24},{36,24},{36,16.98},{19.52,16.98}},color={28,108,
                200}));
        connect(inFlow_liquid.A, hx1simple1.InFlow_liquid) annotation (Line(
              points={{-62.6,19.8},{-62.6,36.9},{-28,36.9},{-28,16.98}},
              color={28,108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-74,-50},{20,-86}},
                lineColor={28,108,200},
                textString="测试换热器元件")}));
      end TestHx1Simple;

      model TestHxCounter
        Components.source_sink_ports.InFlow inFlow(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          m0=-0.5,
          UseT=true,
          p0=1200000,
          T=363.15)
          annotation (Placement(transformation(extent={{-72,-62},{-52,-42}})));
        Components.source_sink_ports.OutFlow outFlow(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          m0=0.5,
          p0=1200000)
          annotation (Placement(transformation(extent={{52,-60},{72,-40}})));
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          N=10,
          Ul=1000,
          Uv=700,
          M=700,
          Cp=4.2e3,
          M_wall=50,
          Cp_wall=0.39e3,
          max_drhode=200,
          witdth_x=0.1,
          Utp=5000,
          Tin_start=363.15,
          Tout_start=388.15,
          T1_liquid=403.15,
          Tn_liquid=373.15)
          annotation (Placement(transformation(extent={{-8,-12},{12,8}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid(T0=403.15)
          annotation (Placement(transformation(extent={{50,10},{70,30}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid(m0=0.5)
          annotation (Placement(transformation(extent={{-64,16},{-44,36}})));
      equation
        connect(inFlow.node, hx1counter.InFlow_working) annotation (Line(points=
               {{-53.8,-52.2},{-53.8,-29.1},{-8,-29.1},{-8,-5.4}}, color={0,0,
                255}));
        connect(hx1counter.OutFlow_working, outFlow.node) annotation (Line(
              points={{11.8,-5.2},{40.9,-5.2},{40.9,-50.2},{70.2,-50.2}}, color=
               {0,0,255}));
        connect(inFlow_liquid.A, hx1counter.InFlow_liquid) annotation (Line(
              points={{69.4,19.8},{69.4,26.9},{11.8,26.9},{11.8,2.6}}, color={
                28,108,200}));
        connect(outFlow_liquid.A, hx1counter.OutFlow_liquid) annotation (Line(
              points={{-64,26},{-36,26},{-36,2.4},{-8,2.4}}, color={28,108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestHxCounter;

      model Test_pumpSimple
        Components.main_components.Pump_simple pump_simple(redeclare package
            medium = ThermoCycle.Media.R245fa_CP,
          p0=250000,
          T0=311.15)
          annotation (Placement(transformation(extent={{-34,-30},{32,26}})));
        Components.source_sink_ports.InFlow inFlow(
          m0=-0.5,
          h0=2.5e5,
          UseT=true,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0=300000,
          T=311.15)
          annotation (Placement(transformation(extent={{-74,-28},{-54,-8}})));
        Components.source_sink_ports.OutFlow outFlow(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          m0=0.5,
          p0=1200000)
          annotation (Placement(transformation(extent={{18,34},{38,54}})));
      equation
        connect(inFlow.node, pump_simple.InFlow) annotation (Line(points={{
                -55.8,-18.2},{-34.9,-18.2},{-34.9,-8.72},{-14.86,-8.72}}, color=
               {0,0,255}));
        connect(outFlow.node, pump_simple.OutFlow) annotation (Line(points={{
                36.2,43.8},{36.2,26.9},{9.56,26.9},{9.56,10.88}}, color={0,0,
                255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end Test_pumpSimple;

      model Test_turbine
        Components.source_sink_ports.InFlow inFlow(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          m0=-0.5,
          p0=1200000,
          UseT=true,
          T=388.15)
          annotation (Placement(transformation(extent={{-70,28},{-50,48}})));
        Components.source_sink_ports.OutFlow outFlow(redeclare package medium
            = ThermoCycle.Media.R245fa_CP, p0=300000)
          annotation (Placement(transformation(extent={{36,26},{56,46}})));
        Components.main_components.Turbine_simple turbine_simple(redeclare
            package medium = ThermoCycle.Media.R245fa_CP, V_s=2.716e-4)
          annotation (Placement(transformation(extent={{-22,0},{-2,20}})));
      equation
        connect(inFlow.node, turbine_simple.InFlow) annotation (Line(points={{
                -51.8,37.8},{-34.9,37.8},{-34.9,14.2},{-17,14.2}}, color={0,0,
                255}));
        connect(turbine_simple.OutFlow, outFlow.node) annotation (Line(points={
                {-7.4,3.4},{70.3,3.4},{70.3,35.8},{54.2,35.8}}, color={0,0,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end Test_turbine;

      model Test_pumptype2
        Components.source_sink_ports.InFlow inFlow(
          m0=-0.5,
          UseT=true,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0=300000,
          T=311.15)
          annotation (Placement(transformation(extent={{-86,-20},{-66,0}})));
        Components.source_sink_ports.OutFlow outFlow(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          m0=0.5,
          p0=1200000)
          annotation (Placement(transformation(extent={{-14,38},{6,58}})));
        Components.main_components.Pump_type2 pump_type2_1(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          T0=311.15,
          v_s=1.2801e-5)
          annotation (Placement(transformation(extent={{-64,-16},{2,42}})));
      equation
        connect(inFlow.node, pump_type2_1.InFlow) annotation (Line(points={{-83.6,
                -10.4},{-55.9,-10.4},{-55.9,6.04},{-44.86,6.04}},       color={
                0,0,255}));
        connect(pump_type2_1.OutFlow, outFlow.node) annotation (Line(points={{
                -20.44,26.34},{-8.22,26.34},{-8.22,47.8},{4.2,47.8}}, color={0,
                0,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end Test_pumptype2;

      model TestTank
        Components.main_components.Tank tank(redeclare package medium =
              ThermoCycle.Media.R245fa_CP, T0(displayUnit="degC") = 313.15)
          annotation (Placement(transformation(extent={{-44,-30},{14,30}})));
        Components.source_sink_ports.InFlow inFlow(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          m0=-1,
          p0=300000,
          UseT=true,
          T=313.15)
          annotation (Placement(transformation(extent={{48,24},{68,44}})));
        Components.source_sink_ports.InFlow inFlow1(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          UseT=true,
          m0=1.05,
          p0=300000,
          T=312.15)
          annotation (Placement(transformation(extent={{-96,22},{-76,42}})));
        Components.units.CellChoiceMedia cellChoiceMedia
          annotation (Placement(transformation(extent={{-116,-122},{-40,-42}})));
      equation
        connect(inFlow.node, tank.InFlow) annotation (Line(points={{66.2,33.8},
                {66.2,11.9},{2.4,11.9},{2.4,12}}, color={0,0,255}));
        connect(inFlow1.node, tank.OutFlow) annotation (Line(points={{-77.8,
                31.8},{-42.9,31.8},{-42.9,-15},{-30.66,-15}}, color={0,0,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestTank;

      model TestCollector
        Components.main_components.SolarCollector solarCollector(
          N=10,
          Area=15,
          K=1,
          M=1,
          Cp=2.1e3,
          m_flow=0.5,
          T0=343.15)
          annotation (Placement(transformation(extent={{-42,-6},{46,56}})));
        Modelica.Blocks.Sources.Constant const(k=700)
          annotation (Placement(transformation(extent={{-44,74},{-34,84}})));
        Modelica.Blocks.Sources.Constant const1(k=25)
          annotation (Placement(transformation(extent={{2,76},{12,86}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid(T0=343.15)
          annotation (Placement(transformation(extent={{-82,-66},{-62,-46}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid(m0=0.5)
          annotation (Placement(transformation(extent={{88,-70},{108,-50}})));
        Components.main_components.Pump_liquid pump_liquid(m0=0.5)
          annotation (Placement(transformation(extent={{60,-8},{80,12}})));
      equation
        connect(const.y, solarCollector.I_in) annotation (Line(points={{-33.5,
                79},{-20,79},{-20,42.98},{-3.28,42.98}}, color={0,0,127}));
        connect(const1.y, solarCollector.Tam_in) annotation (Line(points={{12.5,
                81},{22,81},{22,42.98},{22.24,42.98}}, color={0,0,127}));
        connect(solarCollector.outFlow, pump_liquid.inFlow) annotation (Line(
              points={{46,25.62},{46,23.81},{60,23.81},{60,2}}, color={28,108,
                200}));
        connect(inFlow_liquid.A, outFlow_liquid.A) annotation (Line(points={{
                -62.6,-56.2},{13.7,-56.2},{13.7,-60},{88,-60}}, color={28,108,
                200}));
        connect(pump_liquid.outFlow, solarCollector.inFlow) annotation (Line(
              points={{80,2},{80,-18},{-62,-18},{-62,25},{-42,25}}, color={28,
                108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestCollector;

      model TestValve
        Components.source_sink_ports.InFlow_Tmflow inFlow_Tmflow(
          use_m0=true,
          m0=-1,
          T0=343.15)
          annotation (Placement(transformation(extent={{-90,-10},{-70,10}})));
        Components.source_sink_ports.OutFlow_Tmflow outFlow_Tmflow(use_m0=false,
            m0=1)
          annotation (Placement(transformation(extent={{62,6},{82,26}})));
        Components.source_sink_ports.OutFlow_Tmflow outFlow_Tmflow1(use_m0=
              false, m0=1)
          annotation (Placement(transformation(extent={{46,-52},{66,-32}})));
        Modelica.Blocks.Sources.Constant const1(k=0.6)
          annotation (Placement(transformation(extent={{-26,48},{-16,58}})));
        Components.main_components.Valve valve
          annotation (Placement(transformation(extent={{-14,-8},{6,12}})));
      equation
        connect(inFlow_Tmflow.outFlow, valve.inFlow) annotation (Line(points={{
                -70,0},{-42,0},{-42,2},{-14,2}}, color={28,108,200}));
        connect(const1.y, valve.u) annotation (Line(points={{-15.5,53},{-15.5,
                31.5},{-4,31.5},{-4,7.6}}, color={0,0,127}));
        connect(valve.outFlow_1, outFlow_Tmflow.inFlow) annotation (Line(points=
               {{6,2.2},{34,2.2},{34,15.8},{62,15.8}}, color={28,108,200}));
        connect(valve.outFlow_2, outFlow_Tmflow1.inFlow) annotation (Line(
              points={{-4,-5},{22,-5},{22,-42.2},{46,-42.2}}, color={28,108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestValve;

      model TestPump_liquid
        Components.main_components.Pump_liquid pump_liquid
          annotation (Placement(transformation(extent={{-14,-8},{6,12}})));
        Components.source_sink_ports.InFlow_Tmflow inFlow_Tmflow(
          m0=-0.5,
          T0=363.15,
          use_m0=false)
          annotation (Placement(transformation(extent={{-82,-8},{-62,12}})));
        Components.source_sink_ports.OutFlow_Tmflow outFlow_Tmflow(m0=2, use_m0=
             false)
          annotation (Placement(transformation(extent={{46,-8},{66,12}})));
        Modelica.Blocks.Sources.Constant const(k=0.5)
          annotation (Placement(transformation(extent={{-58,34},{-38,54}})));
      equation
        connect(inFlow_Tmflow.outFlow, pump_liquid.inFlow) annotation (Line(
              points={{-62,2},{-38,2},{-14,2}}, color={28,108,200}));
        connect(pump_liquid.outFlow, outFlow_Tmflow.inFlow) annotation (Line(
              points={{6,2},{46,2},{46,1.8}}, color={28,108,200}));
        connect(const.y, pump_liquid.u) annotation (Line(points={{-37,44},{-22,
                44},{-22,5.9},{-3.9,5.9}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestPump_liquid;

      model TestTank_pump_liquid
        Components.main_components.Tank_liquid tank_liquid(
          Area=1,
          rou=1e3,
          Cp=2.1e3,
          T0=343.15)
          annotation (Placement(transformation(extent={{38,-18},{58,2}})));
        Components.main_components.Pump_liquid pump_liquid
          annotation (Placement(transformation(extent={{24,-42},{4,-22}})));
        Components.source_sink_ports.OutFlow_Tmflow outFlow_Tmflow(m0=1, use_m0=
             false)
          annotation (Placement(transformation(extent={{-40,-44},{-60,-24}})));
        Components.source_sink_ports.InFlow_Tmflow inFlow_Tmflow(
          m0=-0.5,
          use_m0=true,
          T0=343.15)
          annotation (Placement(transformation(extent={{20,28},{40,48}})));
        Modelica.Blocks.Sources.Constant const(k=0.89)
          annotation (Placement(transformation(extent={{-32,-6},{-12,14}})));
      equation
        connect(tank_liquid.outFlow, pump_liquid.inFlow) annotation (Line(
              points={{43.2,-13.6},{33.6,-13.6},{33.6,-32},{24,-32}}, color={28,
                108,200}));
        connect(outFlow_Tmflow.inFlow, pump_liquid.outFlow) annotation (Line(
              points={{-40,-34.2},{-18,-34.2},{-18,-32},{4,-32}}, color={28,108,
                200}));
        connect(inFlow_Tmflow.outFlow, tank_liquid.inFlow) annotation (Line(
              points={{40,38},{68,38},{68,-3.8},{54.2,-3.8}}, color={28,108,200}));
        connect(const.y, pump_liquid.u) annotation (Line(points={{-11,4},{2,4},
                {2,-28.1},{13.9,-28.1}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestTank_pump_liquid;

      model TestMixer
        Components.main_components.mixer mixer
          annotation (Placement(transformation(extent={{-28,-6},{-8,14}})));
        Components.source_sink_ports.InFlow_Tmflow inFlow_Tmflow(
          use_m0=true,
          m0=-1,
          T0=353.15)
          annotation (Placement(transformation(extent={{-80,0},{-60,20}})));
        Components.source_sink_ports.OutFlow_Tmflow outFlow_Tmflow(use_m0=false,
            m0=1)
          annotation (Placement(transformation(extent={{72,16},{92,36}})));
        Components.source_sink_ports.InFlow_Tmflow inFlow_Tmflow1(
          use_m0=true,
          m0=-1,
          T0=343.15)
          annotation (Placement(transformation(extent={{-34,34},{-14,54}})));
      equation
        connect(inFlow_Tmflow1.outFlow, mixer.inFlow2) annotation (Line(points=
                {{-14,44},{-4,44},{-4,7.4},{-12.2,7.4}}, color={28,108,200}));
        connect(inFlow_Tmflow.outFlow, mixer.inFlow1) annotation (Line(points={
                {-60,10},{-42,10},{-42,7.4},{-23.2,7.4}}, color={28,108,200}));
        connect(mixer.outFlow, outFlow_Tmflow.inFlow) annotation (Line(points={
                {-18,1.4},{28,1.4},{28,25.8},{72,25.8}}, color={28,108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestMixer;

      model TestController1
        Components.Controls.PIDcontroller pIDcontroller(u0=0.01, time0=50)
          annotation (Placement(transformation(extent={{-26,-30},{42,60}})));
        Modelica.Blocks.Sources.Constant const4(k=150 + 273)
          annotation (Placement(transformation(extent={{-58,30},{-68,40}})));
        Modelica.Blocks.Sources.Constant const1(k=432)
          annotation (Placement(transformation(extent={{-44,60},{-54,70}})));
      equation
        connect(const4.y, pIDcontroller.u) annotation (Line(points={{-68.5,35},
                {-54.25,35},{-54.25,5.1},{-14.44,5.1}}, color={0,0,127}));
        connect(const1.y, pIDcontroller.u1) annotation (Line(points={{-54.5,65},
                {-37.25,65},{-37.25,28.5},{-14.44,28.5}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestController1;

      model TestDimix
        parameter Integer N=4;
        TJUthermo.Interface.DeMix_node demix(final N=N);
        TJUthermo.Interface.DeMix_node mix(final N=N);
        Interface.SimpleFlange_a In
          annotation (Placement(transformation(extent={{-88,-8},{-68,12}})));
        Interface.SimpleFlange_b Out
          annotation (Placement(transformation(extent={{76,-6},{96,14}})));
      equation
        connect(In,demix.single);
        connect(demix.multi,mix.multi);
        connect(mix.single,Out);

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestDimix;

      model TestTestDimix
        TestDimix testDimix(N=1)
          annotation (Placement(transformation(extent={{-28,-8},{26,40}})));
        Components.source_sink_ports.InFlow_Tmflow inFlow_Tmflow(
          T0=293.15,
          m0=-10,
          use_m0=true)
          annotation (Placement(transformation(extent={{-86,6},{-66,26}})));
        Components.source_sink_ports.OutFlow_Tmflow outFlow_Tmflow(m0=1, use_m0=
             false)
          annotation (Placement(transformation(extent={{70,4},{90,24}})));
      equation
        connect(inFlow_Tmflow.outFlow, testDimix.In) annotation (Line(points={{
                -66,16},{-36,16},{-36,16.48},{-22.06,16.48}}, color={28,108,200}));
        connect(testDimix.Out, outFlow_Tmflow.inFlow) annotation (Line(points={
                {22.22,16.96},{54.11,16.96},{54.11,13.8},{70,13.8}}, color={28,
                108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestTestDimix;

      model TestCrossAir
        Components.units.CrossLiquid crossLiquid(
          F=2,
          V=1,
          U=500,
          Cp=1000,
          M=2,
          N=10)
          annotation (Placement(transformation(extent={{-48,90},{58,12}})));
        Components.units.veryFlow veryFlow(
          N=10,
          L=15,
          dh=0.03,
          m_flow_start=5,
          Tin_start=328.15,
          Tout_start=303.15,
          pstart=300000)
          annotation (Placement(transformation(extent={{-52,-52},{60,-6}})));
        Components.source_sink_ports.DoubleWalls doubleWalls(
          N=10,
          T1=293.15,
          Tn=323.15,
          M=3,
          Cp=500)
          annotation (Placement(transformation(extent={{-36,-10},{56,26}})));
        Components.source_sink_ports.InFlow inFlow(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          m0=-5,
          UseT=true,
          p0(displayUnit="kPa") = 600000,
          T=328.15)
          annotation (Placement(transformation(extent={{-92,-40},{-72,-20}})));
        Components.source_sink_ports.OutFlow outFlow(redeclare package medium
            = ThermoCycle.Media.R245fa_CP, p0=600000)
          annotation (Placement(transformation(extent={{46,-78},{66,-58}})));
        Components.source_sink_ports.MultiSource_liquid multiSource_liquid(
          N=10,
          use_m0=true,
          T0=293.15,
          m0=-10)
          annotation (Placement(transformation(extent={{-96,40},{-76,60}})));
        Components.source_sink_ports.MultiSource_liquid_sink
          multiSource_liquid_sink(
          N=10,
          m0=20,
          use_m0=false)
          annotation (Placement(transformation(extent={{76,42},{96,62}})));
      equation
        connect(doubleWalls.Wall_side1, crossLiquid.wall) annotation (Line(
              points={{10,11.96},{18,11.96},{18,42.81},{4.47,42.81}}, color={28,
                108,200}));
        connect(doubleWalls.Wall_side2, veryFlow.Wall) annotation (Line(points=
                {{10,7.1},{18,7.1},{18,-17.96},{5.12,-17.96}}, color={28,108,
                200}));
        connect(inFlow.node, veryFlow.InFlow) annotation (Line(points={{-73.8,
                -30.2},{-61.9,-30.2},{-61.9,-28.08},{-49.76,-28.08}}, color={0,
                0,255}));
        connect(veryFlow.OutFlow, outFlow.node) annotation (Line(points={{58.88,
                -27.62},{58.88,-45.81},{64.2,-45.81},{64.2,-68.2}}, color={0,0,
                255}));
        connect(multiSource_liquid.outFlow, crossLiquid.A) annotation (Line(
              points={{-76,50},{-46.94,50},{-46.94,51}}, color={28,108,200}));
        connect(crossLiquid.B, multiSource_liquid_sink.inFlow) annotation (Line(
              points={{56.94,51},{67.47,51},{67.47,52},{76.4,52}}, color={28,
                108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestCrossAir;

      model TestCF_liquid
        Components.units.CrossFlow_liquid crossFlow_liquid(
          N=10,
          F=1,
          V=1,
          U=300,
          Cp=1000,
          M=20)
          annotation (Placement(transformation(extent={{-58,-56},{54,0}})));
        Components.units.veryFlow veryFlow(
          N=10,
          L=10,
          dh=0.03,
          Tout_start=323.15,
          pstart=500000)
          annotation (Placement(transformation(extent={{-56,78},{42,18}})));
        Components.source_sink_ports.DoubleWalls doubleWalls(
          N=10,
          T1=293.15,
          Tn=303.15,
          M=20,
          Cp=1000)
          annotation (Placement(transformation(extent={{-42,-12},{38,24}})));
        Components.source_sink_ports.MultiSource_liquid multiSource_liquid(
          m0=-2,
          use_m0=true,
          N=10,
          T0=293.15)
          annotation (Placement(transformation(extent={{-106,-40},{-86,-20}})));
        Components.source_sink_ports.MultiSource_liquid_sink
          multiSource_liquid_sink(
          N=10,
          m0=2,
          use_m0=false)
          annotation (Placement(transformation(extent={{70,-38},{90,-18}})));
        Components.source_sink_ports.InFlow inFlow(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          UseT=true,
          p0(displayUnit="kPa") = 600000,
          m0=-2,
          T=328.15)
          annotation (Placement(transformation(extent={{-98,38},{-78,58}})));
        Components.source_sink_ports.OutFlow outFlow1(
                                                     redeclare package medium
            = ThermoCycle.Media.R245fa_CP, p0=600000)
          annotation (Placement(transformation(extent={{88,40},{68,60}})));
      equation
        connect(doubleWalls.Wall_side1, veryFlow.Wall) annotation (Line(points=
                {{-2,9.96},{-6,9.96},{-6,30},{-6.02,30},{-6.02,33.6}}, color={
                28,108,200}));
        connect(doubleWalls.Wall_side2, crossFlow_liquid.node_wall) annotation (
           Line(points={{-2,5.1},{-3.68,5.1},{-3.68,-19.32}}, color={28,108,200}));
        connect(multiSource_liquid.outFlow, crossFlow_liquid.node_a)
          annotation (Line(points={{-86,-30},{-74,-30},{-74,-28},{-58,-28}},
              color={28,108,200}));
        connect(crossFlow_liquid.node_b, multiSource_liquid_sink.inFlow)
          annotation (Line(points={{52.88,-27.44},{62.44,-27.44},{62.44,-28},{
                70.4,-28}}, color={28,108,200}));
        connect(inFlow.node, veryFlow.InFlow) annotation (Line(points={{-79.8,
                47.8},{-66.9,47.8},{-66.9,46.8},{-54.04,46.8}}, color={0,0,255}));
        connect(veryFlow.OutFlow, outFlow1.node) annotation (Line(points={{
                41.02,46.2},{55.51,46.2},{55.51,49.8},{69.8,49.8}}, color={0,0,
                255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestCF_liquid;

      model TestHXCross
        Components.main_components.Hx1CrossFlow hx1CrossFlow(
          N=10,
          M_wall=20,
          Cp_wall=500,
          np=8,
          U=172.6,
          F_liquid=0.1396,
          L=0.365,
          dh=0.0015,
          Cp=1017,
          m_working=0.05,
          m_liquid=0.2816,
          pstart(displayUnit="kPa") = 1682000,
          V_liquid=0.00307,
          M=0.0031,
          Tin_start=347.15,
          Tout_start=347.15,
          T1_liquid=298.15,
          Tn_liquid=318.15)
          annotation (Placement(transformation(extent={{-54,-38},{46,42}})));
        Components.source_sink_ports.InFlow inFlow(
          UseT=true,
          redeclare package medium = ThermoCycle.Media.R134a_CP,
          m0=-0.05,
          p0(displayUnit="kPa") = 1682000,
          T(displayUnit="K") = 358)
          annotation (Placement(transformation(extent={{-96,-10},{-76,10}})));
        Components.source_sink_ports.OutFlow outFlow(redeclare package medium
            = ThermoCycle.Media.R134a_CP, p0(displayUnit="kPa") = 1682000)
          annotation (Placement(transformation(extent={{88,-16},{68,4}})));
        Components.source_sink_ports.MultiSource_liquid multiSource_liquid(
          N=10,
          use_m0=true,
          m0=-0.2816,
          T0=298.15)
          annotation (Placement(transformation(extent={{-40,-66},{-20,-46}})));
        Components.source_sink_ports.MultiSource_liquid_sink
          multiSource_liquid_sink(
          N=10,
          m0=20,
          use_m0=false)
          annotation (Placement(transformation(extent={{16,50},{36,70}})));
      equation
        connect(inFlow.node, hx1CrossFlow.InFlow) annotation (Line(points={{
                -77.8,-0.2},{-65.9,-0.2},{-65.9,0.4},{-55,0.4}}, color={0,0,255}));
        connect(hx1CrossFlow.OutFlow, outFlow.node) annotation (Line(points={{
                47,1.2},{60.5,1.2},{60.5,-6.2},{69.8,-6.2}}, color={0,0,255}));
        connect(multiSource_liquid.outFlow, hx1CrossFlow.LiquidIn) annotation (
            Line(points={{-20,-56},{-12,-56},{-12,-30.8},{-4,-30.8}}, color={28,
                108,200}));
        connect(multiSource_liquid_sink.inFlow, hx1CrossFlow.LiquidOut)
          annotation (Line(points={{16.4,60},{6,60},{6,32.4},{-4,32.4}}, color=
                {28,108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestHXCross;

      model Test_new_collectors
        Modelica.Blocks.Sources.Constant const(k=700)
          annotation (Placement(transformation(extent={{-34,84},{-24,94}})));
        Modelica.Blocks.Sources.Constant const1(k=25)
          annotation (Placement(transformation(extent={{12,86},{22,96}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid(T0=343.15)
          annotation (Placement(transformation(extent={{-72,-56},{-52,-36}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid(m0=1)
          annotation (Placement(transformation(extent={{98,-60},{118,-40}})));
        Components.main_components.SolarCollector_New solarCollector_New(
          N=10,
          Area=15,
          K=1,
          M=10,
          Cp=2400,
          m_flow=1,
          ti=10,
          lati=35,
          Date=100,
          T0=343.15,
          use_eff2=false)
          annotation (Placement(transformation(extent={{-46,10},{42,72}})));
      equation
        connect(const.y, solarCollector_New.I_in) annotation (Line(points={{
                -23.5,89},{-23.5,73.5},{-7.28,73.5},{-7.28,58.98}}, color={0,0,
                127}));
        connect(const1.y, solarCollector_New.Tam_in) annotation (Line(points={{
                22.5,91},{22.5,74.5},{18.24,74.5},{18.24,58.98}}, color={0,0,
                127}));
        connect(inFlow_liquid.A, solarCollector_New.inFlow) annotation (Line(
              points={{-52.6,-46.2},{-52.6,-3.1},{-46,-3.1},{-46,41}}, color={
                28,108,200}));
        connect(solarCollector_New.outFlow, outFlow_liquid.A) annotation (Line(
              points={{42,41.62},{72,41.62},{72,-50},{98,-50}}, color={28,108,
                200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end Test_new_collectors;
    end TestComponets;

    package TestCycle

      model TestSimpleORC_1
        Components.source_sink_ports.InFlow_liquid inFlow_liquid(T0=403.15)
          annotation (Placement(transformation(extent={{62,40},{82,60}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid(m0=0.5)
          annotation (Placement(transformation(extent={{-58,38},{-38,58}})));
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Uv=700,
          Cp=4.2e3,
          Cp_wall=0.39e3,
          max_drhode=200,
          witdth_x=0.1,
          Utp=5000,
          F_working=8.8,
          F_liquid=8.8,
          N=20,
          M_wall=10,
          M=300,
          Tin_start=312.15,
          Tout_start=388.15,
          T1_liquid=403.15,
          Tn_liquid=373.15)
          annotation (Placement(transformation(extent={{-14,22},{6,42}})));
        Components.main_components.Turbine_simple turbine_simple(redeclare
            package medium = ThermoCycle.Media.R245fa_CP,
          V_s=2.716e-4,
          Tin_start=383.15)
          annotation (Placement(transformation(extent={{30,-16},{50,4}})));
        Components.main_components.Pump_simple pump_simple(redeclare package
            medium = ThermoCycle.Media.R245fa_CP,
          p0=300000,
          T0=311.15)
          annotation (Placement(transformation(extent={{-82,-34},{-16,22}})));
        Components.main_components.Hx1counter condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          N=10,
          Ul=1000,
          Uv=700,
          M=700,
          Cp=4.2e3,
          M_wall=50,
          Cp_wall=0.39e3,
          witdth_x=0.1,
          Utp=5000,
          F_working=13.8,
          F_liquid=13.8,
          m_working=0.5,
          max_drhode=98,
          Tin_start=371.15,
          Tout_start=303.15,
          T1_liquid=293.15,
          Tn_liquid=310.15,
          m_liquid=3.6,
          pstart=300000)
          annotation (Placement(transformation(extent={{22,-46},{-26,-82}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid1(T0=293.15)
          annotation (Placement(transformation(extent={{-66,-104},{-46,-84}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(m0=3.5)
          annotation (Placement(transformation(extent={{56,-98},{76,-78}})));
        Components.source_sink_ports.InFlow inFlow(
          m0=-0.5,
          UseT=true,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0=300000,
          T=311.15)
          annotation (Placement(transformation(extent={{-96,-30},{-76,-10}})));
        Components.source_sink_ports.OutFlow outFlow(m0=0.5,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0=300000)
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=0,
              origin={-80,-60})));
      equation
        connect(outFlow_liquid.A, hx1counter.OutFlow_liquid) annotation (Line(
              points={{-58,48},{-28,48},{-28,36.4},{-14,36.4}}, color={28,108,
                200}));
        connect(inFlow_liquid.A, hx1counter.InFlow_liquid) annotation (Line(
              points={{81.4,49.8},{51.7,49.8},{51.7,36.6},{5.8,36.6}}, color={
                28,108,200}));
        connect(hx1counter.OutFlow_working, turbine_simple.InFlow) annotation (
            Line(points={{5.8,28.8},{20.9,28.8},{20.9,-1.8},{35,-1.8}}, color={
                0,0,255}));
        connect(pump_simple.OutFlow, hx1counter.InFlow_working) annotation (
            Line(points={{-38.44,6.88},{-43.22,6.88},{-43.22,28.6},{-14,28.6}},
                        color={0,0,255}));
        connect(inFlow_liquid1.A, condenser.InFlow_liquid) annotation (Line(
              points={{-46.6,-94.2},{-33.3,-94.2},{-33.3,-72.28},{-25.52,-72.28}},
              color={28,108,200}));
        connect(condenser.OutFlow_liquid, outFlow_liquid1.A) annotation (Line(
              points={{22,-71.92},{40,-71.92},{40,-88},{56,-88}}, color={28,108,
                200}));
        connect(turbine_simple.OutFlow, condenser.InFlow_working) annotation (
            Line(points={{44.6,-12.6},{44.6,-38.3},{22,-38.3},{22,-57.88}},
              color={0,0,255}));
        connect(inFlow.node, pump_simple.InFlow) annotation (Line(points={{-93.6,
                -20.4},{-71.9,-20.4},{-71.9,-12.72},{-62.86,-12.72}},
              color={0,0,255}));
        connect(outFlow.node, condenser.OutFlow_working) annotation (Line(
              points={{-88.2,-60.2},{-57.1,-60.2},{-57.1,-58.24},{-25.52,-58.24}},
              color={0,0,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-40,88},{46,64}},
                lineColor={28,108,200},
                textString="initial problem (暂未解决)")}));
      end TestSimpleORC_1;

      model TestSimpleORC_2
        Components.source_sink_ports.InFlow_liquid inFlow_liquid(T0=403.15)
          annotation (Placement(transformation(extent={{66,56},{86,76}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid(m0=0.5)
          annotation (Placement(transformation(extent={{-54,54},{-34,74}})));
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Uv=700,
          Cp=4.2e3,
          Cp_wall=0.39e3,
          max_drhode=200,
          witdth_x=0.1,
          Utp=5000,
          F_working=8.8,
          F_liquid=8.8,
          N=20,
          M_wall=10,
          M=300,
          Tin_start=312.15,
          Tout_start=388.15,
          T1_liquid=403.15,
          Tn_liquid=373.15)
          annotation (Placement(transformation(extent={{-18,30},{24,64}})));
        Components.main_components.Turbine_simple turbine_simple(redeclare
            package medium = ThermoCycle.Media.R245fa_CP,
          V_s=2.716e-4,
          Tin_start=383.15)
          annotation (Placement(transformation(extent={{34,0},{54,20}})));
        Components.main_components.Hx1counter condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          N=10,
          Ul=1000,
          Uv=700,
          M=700,
          Cp=4.2e3,
          M_wall=50,
          Cp_wall=0.39e3,
          witdth_x=0.1,
          Utp=5000,
          F_working=13.8,
          F_liquid=13.8,
          m_working=0.5,
          m_liquid=3.6,
          max_drhode=200,
          Tin_start=371.15,
          Tout_start=303.15,
          T1_liquid=303.15,
          Tn_liquid=313.15,
          pstart=300000)
          annotation (Placement(transformation(extent={{26,-30},{-22,-66}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid1(T0=303.15)
          annotation (Placement(transformation(extent={{-46,-92},{-26,-72}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(m0=3.5)
          annotation (Placement(transformation(extent={{60,-82},{80,-62}})));
        Components.source_sink_ports.InFlow inFlow(
          m0=-0.5,
          UseT=true,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0=300000,
          T=311.15)
          annotation (Placement(transformation(extent={{-92,-14},{-72,6}})));
        Components.source_sink_ports.OutFlow outFlow(m0=0.5,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0=300000)
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=0,
              origin={-76,-44})));
        Components.main_components.Pump_simple pump_simple1(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          steadytime=5,
          T0=311.15,
          steadyInFlow=false)
          annotation (Placement(transformation(extent={{-68,-14},{-2,38}})));
      equation
        connect(outFlow_liquid.A,hx1counter. OutFlow_liquid) annotation (Line(
              points={{-54,64},{-24,64},{-24,54.48},{-18,54.48}},
                                                                color={28,108,
                200}));
        connect(inFlow_liquid.A,hx1counter. InFlow_liquid) annotation (Line(
              points={{85.4,65.8},{55.7,65.8},{55.7,54.82},{23.58,54.82}},
                                                                       color={
                28,108,200}));
        connect(hx1counter.OutFlow_working,turbine_simple. InFlow) annotation (
            Line(points={{23.58,41.56},{24.9,41.56},{24.9,14.2},{39,14.2}},
                                                                        color={
                0,0,255}));
        connect(inFlow_liquid1.A,condenser. InFlow_liquid) annotation (Line(
              points={{-26.6,-82.2},{-29.3,-82.2},{-29.3,-56.28},{-21.52,-56.28}},
              color={28,108,200}));
        connect(condenser.OutFlow_liquid,outFlow_liquid1. A) annotation (Line(
              points={{26,-55.92},{44,-55.92},{44,-72},{60,-72}}, color={28,108,
                200}));
        connect(turbine_simple.OutFlow,condenser. InFlow_working) annotation (
            Line(points={{48.6,3.4},{48.6,-22.3},{26,-22.3},{26,-41.88}},
              color={0,0,255}));
        connect(inFlow.node, outFlow.node) annotation (Line(points={{-89.6,-4.4},
                {-89.6,-24.1},{-84.2,-24.1},{-84.2,-44.2}}, color={0,0,255}));
        connect(pump_simple1.InFlow, condenser.OutFlow_working) annotation (
            Line(points={{-48.86,5.76},{-48.86,-19.12},{-21.52,-19.12},{-21.52,
                -42.24}}, color={0,0,255}));
        connect(pump_simple1.OutFlow, hx1counter.InFlow_working) annotation (
            Line(points={{-24.44,23.96},{-24.44,33.98},{-18,33.98},{-18,41.22}},
              color={0,0,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-36,104},{50,80}},
                lineColor={28,108,200},
                textString="initial problem (已解决  强硬设置初值)")}));
      end TestSimpleORC_2;

      model TestSimpleORCwithType2
        Components.source_sink_ports.InFlow_liquid inFlow_liquid(T0=403.15)
          annotation (Placement(transformation(extent={{72,50},{92,70}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid(m0=0.5)
          annotation (Placement(transformation(extent={{-48,48},{-28,68}})));
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Uv=700,
          Cp=4.2e3,
          Cp_wall=0.39e3,
          max_drhode=200,
          witdth_x=0.1,
          Utp=5000,
          F_working=8.8,
          F_liquid=8.8,
          N=20,
          M_wall=10,
          M=300,
          Tin_start=312.15,
          Tout_start=388.15,
          T1_liquid=403.15,
          Tn_liquid=373.15)
          annotation (Placement(transformation(extent={{-4,32},{16,52}})));
        Components.main_components.Turbine_simple turbine_simple(redeclare
            package medium = ThermoCycle.Media.R245fa_CP,
          V_s=2.716e-4,
          Tin_start=383.15)
          annotation (Placement(transformation(extent={{40,-6},{60,14}})));
        Components.main_components.Hx1counter condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          N=10,
          Ul=1000,
          Uv=700,
          M=700,
          Cp=4.2e3,
          M_wall=50,
          Cp_wall=0.39e3,
          witdth_x=0.1,
          Utp=5000,
          F_working=13.8,
          F_liquid=13.8,
          m_working=0.5,
          max_drhode=98,
          Tin_start=371.15,
          Tout_start=303.15,
          T1_liquid=293.15,
          Tn_liquid=310.15,
          m_liquid=3.6,
          pstart=300000)
          annotation (Placement(transformation(extent={{32,-36},{-16,-72}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid1(T0=293.15)
          annotation (Placement(transformation(extent={{-56,-94},{-36,-74}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(m0=3.5)
          annotation (Placement(transformation(extent={{66,-88},{86,-68}})));
        Components.main_components.Pump_type2 pump_type2_1(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          v_s=1.2801e-5,
          T0=311.15)
          annotation (Placement(transformation(extent={{-70,-16},{-4,42}})));
      equation
        connect(outFlow_liquid.A,hx1counter. OutFlow_liquid) annotation (Line(
              points={{-48,58},{-18,58},{-18,46.4},{-4,46.4}},  color={28,108,
                200}));
        connect(inFlow_liquid.A,hx1counter. InFlow_liquid) annotation (Line(
              points={{91.4,59.8},{61.7,59.8},{61.7,46.6},{15.8,46.6}},color={
                28,108,200}));
        connect(hx1counter.OutFlow_working,turbine_simple. InFlow) annotation (
            Line(points={{15.8,38.8},{30.9,38.8},{30.9,8.2},{45,8.2}},  color={
                0,0,255}));
        connect(inFlow_liquid1.A,condenser. InFlow_liquid) annotation (Line(
              points={{-36.6,-84.2},{-23.3,-84.2},{-23.3,-62.28},{-15.52,-62.28}},
              color={28,108,200}));
        connect(condenser.OutFlow_liquid,outFlow_liquid1. A) annotation (Line(
              points={{32,-61.92},{50,-61.92},{50,-78},{66,-78}}, color={28,108,
                200}));
        connect(turbine_simple.OutFlow,condenser. InFlow_working) annotation (
            Line(points={{54.6,-2.6},{54.6,-28.3},{32,-28.3},{32,-47.88}},
              color={0,0,255}));
        connect(pump_type2_1.OutFlow, hx1counter.InFlow_working) annotation (
            Line(points={{-26.44,26.34},{-15.22,26.34},{-15.22,38.6},{-4,38.6}},
              color={0,0,255}));
        connect(pump_type2_1.InFlow, condenser.OutFlow_working) annotation (
            Line(points={{-50.86,6.04},{-50.86,-20.98},{-15.52,-20.98},{-15.52,
                -48.24}}, color={0,0,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-30,98},{56,74}},
                lineColor={28,108,200},
                textString="initial problem (暂未解决)")}));
      end TestSimpleORCwithType2;

      model TestSimpleORC_3
        Components.main_components.Tank tank(redeclare package medium =
              ThermoCycle.Media.R245fa_CP, T0(displayUnit="degC") = 313.15,
          Aera=0.5,
          L0=0.5)
          annotation (Placement(transformation(extent={{-84,-74},{-26,-14}})));
        Components.source_sink_ports.InFlow inFlow(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          UseT=true,
          m0=-0.5,
          p0=1200000,
          T=387.15)
          annotation (Placement(transformation(extent={{52,26},{72,46}})));
        Components.source_sink_ports.OutFlow outFlow(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          m0=0.5,
          p0=1200000)
          annotation (Placement(transformation(extent={{62,48},{82,68}})));
        Components.main_components.Hx1counter condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          N=10,
          Ul=1000,
          Uv=700,
          M=700,
          Cp=4.2e3,
          M_wall=50,
          Cp_wall=0.39e3,
          witdth_x=0.1,
          Utp=5000,
          F_working=13.8,
          F_liquid=13.8,
          m_working=0.5,
          m_liquid=3.6,
          max_drhode=200,
          Tin_start=371.15,
          Tout_start=303.15,
          T1_liquid=303.15,
          Tn_liquid=313.15,
          pstart=300000)
          annotation (Placement(transformation(extent={{48,-40},{0,-76}})));
        Components.main_components.Turbine_simple turbine_simple(redeclare
            package medium = ThermoCycle.Media.R245fa_CP,
          V_s=2.716e-4,
          n_turbine=0.7,
          Tin_start=383.15)
          annotation (Placement(transformation(extent={{52,-26},{92,16}})));
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Uv=700,
          Cp=4.2e3,
          Cp_wall=0.39e3,
          max_drhode=200,
          witdth_x=0.1,
          Utp=5000,
          F_working=8.8,
          F_liquid=8.8,
          N=20,
          M_wall=10,
          M=300,
          Tin_start=312.15,
          Tout_start=388.15,
          T1_liquid=403.15,
          Tn_liquid=373.15)
          annotation (Placement(transformation(extent={{-24,40},{18,74}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid(T0=403.15)
          annotation (Placement(transformation(extent={{72,62},{52,82}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid(m0=0.5)
          annotation (Placement(transformation(extent={{-64,70},{-84,90}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid1(T0=303.15)
          annotation (Placement(transformation(extent={{-48,-94},{-28,-74}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(m0=3.6)
          annotation (Placement(transformation(extent={{68,-76},{88,-56}})));
        Components.main_components.pump_type3 pump_type3_1(
          v_s=1.2801e-5,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          T0=311.15)
          annotation (Placement(transformation(extent={{-36,-34},{-102,34}})));
        Modelica.Blocks.Sources.Step step(
          height=5,
          offset=30,
          startTime=500)
          annotation (Placement(transformation(extent={{-128,24},{-108,44}})));
      equation
        connect(tank.InFlow, condenser.OutFlow_working) annotation (Line(points=
               {{-37.6,-32},{-18,-32},{-18,-52.24},{0.48,-52.24}}, color={0,0,
                255}));
        connect(turbine_simple.OutFlow, condenser.InFlow_working) annotation (
            Line(points={{81.2,-18.86},{81.2,-37.3},{48,-37.3},{48,-51.88}},
              color={0,0,255}));
        connect(inFlow.node, outFlow.node) annotation (Line(points={{54.4,35.6},
                {54.4,46.9},{80.2,46.9},{80.2,57.8}}, color={0,0,255}));
        connect(hx1counter.OutFlow_working, turbine_simple.InFlow) annotation (
            Line(points={{17.58,51.56},{17.58,19.78},{62,19.78},{62,3.82}},
              color={0,0,255}));
        connect(outFlow_liquid.A, hx1counter.OutFlow_liquid) annotation (Line(
              points={{-64,80},{-52,80},{-52,64.48},{-24,64.48}}, color={28,108,
                200}));
        connect(inFlow_liquid.A, hx1counter.InFlow_liquid) annotation (Line(
              points={{52.6,71.8},{49.7,71.8},{49.7,64.82},{17.58,64.82}},
              color={28,108,200}));
        connect(inFlow_liquid1.A, condenser.InFlow_liquid) annotation (Line(
              points={{-28.6,-84.2},{-14.3,-84.2},{-14.3,-66.28},{0.48,-66.28}},
              color={28,108,200}));
        connect(condenser.OutFlow_liquid, outFlow_liquid1.A) annotation (Line(
              points={{48,-65.92},{58,-65.92},{58,-66},{68,-66}}, color={28,108,
                200}));
        connect(pump_type3_1.InFlow, tank.OutFlow) annotation (Line(points={{
                -55.14,-8.16},{-55.14,-14.08},{-70.66,-14.08},{-70.66,-59}},
              color={0,0,255}));
        connect(pump_type3_1.OutFlow, hx1counter.InFlow_working) annotation (
            Line(points={{-79.56,15.64},{-80.78,15.64},{-80.78,51.22},{-24,
                51.22}}, color={0,0,255}));
        connect(step.y, pump_type3_1.u) annotation (Line(points={{-107,34},{
                -107,4},{-78.9,4},{-78.9,4.08}},   color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-40,4},{46,-20}},
                lineColor={28,108,200},
                textString="添加Tank原件后整体系统稳定很多
便于调通")}));
      end TestSimpleORC_3;

      model TestOilCycle
        Components.main_components.Pump_liquid pump_liquid(m0=1)
          annotation (Placement(transformation(extent={{30,-78},{10,-58}})));
        Modelica.Blocks.Sources.Constant const(k=1)
          annotation (Placement(transformation(extent={{2,-48},{12,-38}})));
        Modelica.Blocks.Sources.Constant const1(
                                               k=700)
          annotation (Placement(transformation(extent={{-66,26},{-56,36}})));
        Modelica.Blocks.Sources.Constant const2(k=25)
          annotation (Placement(transformation(extent={{-18,32},{-8,42}})));
        Components.main_components.SolarCollector solarCollector1(
          N=10,
          Area=15,
          K=1,
          M=1,
          Cp=2.1e3,
          m_flow=1,
          T0=343.15)
          annotation (Placement(transformation(extent={{-48,-22},{-10,18}})));
        Components.main_components.Tank_liquid tank_liquid(
          Area=0.5,
          rou=1000,
          Cp=2.1e3,
          m_flow=1,
          T0=343.15)
          annotation (Placement(transformation(extent={{58,-38},{78,-18}})));
      equation
        connect(const.y, pump_liquid.u) annotation (Line(points={{12.5,-43},{20,
                -43},{20,-64.1},{19.9,-64.1}}, color={0,0,127}));
        connect(const2.y, solarCollector1.Tam_in) annotation (Line(points={{
                -7.5,37},{-7.5,23.5},{-20.26,23.5},{-20.26,9.6}}, color={0,0,
                127}));
        connect(solarCollector1.inFlow, pump_liquid.outFlow) annotation (Line(
              points={{-48,-2},{-66,-2},{-66,-68},{10,-68}}, color={28,108,200}));
        connect(const1.y, solarCollector1.I_in) annotation (Line(points={{-55.5,
                31},{-30.75,31},{-30.75,9.6},{-31.28,9.6}}, color={0,0,127}));
        connect(solarCollector1.outFlow, tank_liquid.inFlow) annotation (Line(
              points={{-10,-1.6},{32,-1.6},{32,-23.8},{74.2,-23.8}}, color={28,
                108,200}));
        connect(tank_liquid.outFlow, pump_liquid.inFlow) annotation (Line(
              points={{63.2,-33.6},{63.2,-50.8},{30,-50.8},{30,-68}}, color={28,
                108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestOilCycle;

      model TestOil_ORC
        Components.main_components.Pump_liquid pump_liquid
          annotation (Placement(transformation(extent={{48,-34},{28,-14}})));
        Modelica.Blocks.Sources.Constant const(k=0.5)
          annotation (Placement(transformation(extent={{20,-4},{30,6}})));
        Modelica.Blocks.Sources.Constant const1(
                                               k=700)
          annotation (Placement(transformation(extent={{-48,70},{-38,80}})));
        Modelica.Blocks.Sources.Constant const2(k=25)
          annotation (Placement(transformation(extent={{0,76},{10,86}})));
        Components.main_components.SolarCollector solarCollector1(
          K=1,
          m_flow=0.5,
          N=20,
          M=20,
          Area=300,
          Cp=4.2e3,
          T0=403.15)
          annotation (Placement(transformation(extent={{-30,22},{8,62}})));
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Uv=700,
          Cp=4.2e3,
          Cp_wall=0.39e3,
          max_drhode=200,
          witdth_x=0.1,
          Utp=5000,
          F_working=8.8,
          F_liquid=8.8,
          N=20,
          M_wall=10,
          M=300,
          Tin_start=312.15,
          Tout_start=388.15,
          T1_liquid=403.15,
          Tn_liquid=373.15)
          annotation (Placement(transformation(extent={{-38,-40},{4,-6}})));
        Components.source_sink_ports.InFlow inFlow(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          UseT=true,
          m0=-0.5,
          p0=1200000,
          T=312.15)
          annotation (Placement(transformation(extent={{-92,-70},{-72,-50}})));
        Components.source_sink_ports.OutFlow outFlow(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          m0=0.5,
          p0=1200000)
          annotation (Placement(transformation(extent={{-12,-82},{8,-62}})));
        Components.main_components.Tank_liquid tank_liquid(
          m_flow=0.5,
          Area=0.1,
          rou=1000,
          L0=0.3,
          Cp=4.2e3,
          T0=403.15)
          annotation (Placement(transformation(extent={{54,-10},{74,10}})));
      equation
        connect(const.y,pump_liquid. u) annotation (Line(points={{30.5,1},{38,1},
                {38,-20.1},{37.9,-20.1}},      color={0,0,127}));
        connect(const2.y,solarCollector1. Tam_in) annotation (Line(points={{10.5,81},
                {10.5,67.5},{-2.26,67.5},{-2.26,53.6}},           color={0,0,
                127}));
        connect(const1.y, solarCollector1.I_in) annotation (Line(points={{-37.5,
                75},{-12.75,75},{-12.75,53.6},{-13.28,53.6}}, color={0,0,127}));
        connect(hx1counter.InFlow_liquid, pump_liquid.outFlow) annotation (Line(
              points={{3.58,-15.18},{15.79,-15.18},{15.79,-24},{28,-24}}, color=
               {28,108,200}));
        connect(solarCollector1.inFlow, hx1counter.OutFlow_liquid) annotation (
            Line(points={{-30,42},{-64,42},{-64,-15.52},{-38,-15.52}}, color={
                28,108,200}));
        connect(inFlow.node, hx1counter.InFlow_working) annotation (Line(points={{-73.8,
                -60.2},{-55.9,-60.2},{-55.9,-28.78},{-38,-28.78}},        color=
               {0,0,255}));
        connect(hx1counter.OutFlow_working, outFlow.node) annotation (Line(
              points={{3.58,-28.44},{3.58,-60.22},{6.2,-60.22},{6.2,-72.2}},
              color={0,0,255}));
        connect(solarCollector1.outFlow, tank_liquid.inFlow) annotation (Line(
              points={{8,42.4},{40,42.4},{40,4.2},{70.2,4.2}}, color={28,108,
                200}));
        connect(tank_liquid.outFlow, pump_liquid.inFlow) annotation (Line(
              points={{59.2,-5.6},{59.2,-14.8},{48,-14.8},{48,-24}}, color={28,
                108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-50,22},{-6,-16}},
                lineColor={28,108,200},
                textString="流体比热没有规定成一个值(已解决)")}));
      end TestOil_ORC;

      model TestORCwithCollector
        Components.main_components.Pump_liquid pump_liquid
          annotation (Placement(transformation(extent={{60,4},{40,24}})));
        Modelica.Blocks.Sources.Constant const(k=3.5)
          annotation (Placement(transformation(extent={{46,24},{36,34}})));
        Modelica.Blocks.Sources.Constant const1(k=450)
          annotation (Placement(transformation(extent={{-74,72},{-64,82}})));
        Modelica.Blocks.Sources.Constant const2(k=25)
          annotation (Placement(transformation(extent={{10,86},{20,96}})));
        Components.main_components.SolarCollector solarCollector1(
          K=1,
          N=20,
          Area=300,
          M=137,
          m_flow=4.5,
          Cp=2.34e3,
          T0=403.15)
          annotation (Placement(transformation(extent={{-20,32},{18,72}})));
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Uv=700,
          Cp_wall=0.39e3,
          max_drhode=200,
          witdth_x=0.1,
          Utp=5000,
          F_working=8.8,
          F_liquid=8.8,
          N=20,
          M_wall=10,
          M=300,
          Cp=2.34e3,
          m_liquid=4.5,
          Tin_start=312.15,
          Tout_start=388.15,
          T1_liquid=403.15,
          Tn_liquid=373.15,
          m_working=0.3)
          annotation (Placement(transformation(extent={{-30,-2},{12,32}})));
        Components.main_components.Tank tank(redeclare package medium =
              ThermoCycle.Media.R245fa_CP, T0(displayUnit="degC") = 313.15)
          annotation (Placement(transformation(extent={{-80,-72},{-38,-36}})));
        Components.main_components.Pump_type2 pump_type2_1(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          v_s=1.2801e-5,
          T0=311.15,
          Np=18)
          annotation (Placement(transformation(extent={{-110,-36},{-44,22}})));
        Components.main_components.Hx1counter condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          N=10,
          Ul=1000,
          Uv=700,
          M=700,
          Cp=4.2e3,
          M_wall=50,
          Cp_wall=0.39e3,
          witdth_x=0.1,
          Utp=5000,
          F_working=13.8,
          F_liquid=13.8,
          max_drhode=200,
          m_liquid=5,
          m_working=0.3,
          Tin_start=371.15,
          Tout_start=303.15,
          T1_liquid=303.15,
          Tn_liquid=313.15,
          pstart=300000)
          annotation (Placement(transformation(extent={{36,-46},{-12,-82}})));
        Components.main_components.Turbine_simple turbine_simple(redeclare
            package medium = ThermoCycle.Media.R245fa_CP,
          n_turbine=0.7,
          Np=18,
          V_s=3.121e-4,
          Tin_start=383.15)
          annotation (Placement(transformation(extent={{40,-52},{60,-32}})));
        Components.source_sink_ports.OutFlow_Tmflow outFlow_Tmflow(m0=1, use_m0=
             false)
          annotation (Placement(transformation(extent={{70,-96},{90,-76}})));
        Components.main_components.Tank_liquid tank_liquid1(
          rou=825,
          Cp=2.34e3,
          Area=1,
          L0=0.8,
          T0=373.15)
          annotation (Placement(transformation(extent={{22,66},{42,86}})));
        Components.source_sink_ports.InFlow_Tmflow inFlow_Tmflow(
          use_m0=true,
          T0=298.15,
          m0=-5)
          annotation (Placement(transformation(extent={{-92,-92},{-72,-72}})));
        Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(
          tableOnFile=true,
          columns={2},
          tableName="uaua",
          fileName="uaua.mat")
          annotation (Placement(transformation(extent={{-38,82},{-24,96}})));
        Components.main_components.Valve valve annotation (Placement(
              transformation(
              extent={{-6,6},{6,-6}},
              rotation=-90,
              origin={80,44})));
        Modelica.Blocks.Sources.Constant const3(k=0.99)
          annotation (Placement(transformation(extent={{54,34},{64,44}})));
        Components.main_components.mixer mixer annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-56,50})));
        Components.main_components.Pump_liquid pump_liquid1
          annotation (Placement(transformation(extent={{84,62},{64,82}})));
        Modelica.Blocks.Sources.Constant const4(k=130 + 273)
          annotation (Placement(transformation(extent={{24,108},{14,118}})));
        Components.main_components.sens_liquid sens_liquid
          annotation (Placement(transformation(extent={{36,42},{56,62}})));
        Modelica.Blocks.Math.Add add(k1=-1)
          annotation (Placement(transformation(extent={{50,102},{60,112}})));
        Components.Controls.TimeModel timeModel(
          umax=100,
          umin=1,
          dt=2800)
          annotation (Placement(transformation(extent={{70,116},{90,136}})));
        Components.main_components.Tank_liquid tank_liquid2(
          rou=825,
          Cp=2.34e3,
          Area=0.1,
          L0=0.1,
          T0=403.15)
          annotation (Placement(transformation(extent={{68,-14},{88,6}})));
        Modelica.Blocks.Math.Gain gain(k=0.08)
          annotation (Placement(transformation(extent={{118,110},{126,118}})));
        Modelica.Blocks.Sources.Constant const5(k=0.01)
          annotation (Placement(transformation(extent={{104,88},{94,98}})));
      equation
        connect(const.y,pump_liquid. u) annotation (Line(points={{35.5,29},{
                34.5,29},{34.5,17.9},{49.9,17.9}},
                                               color={0,0,127}));
        connect(const2.y,solarCollector1. Tam_in) annotation (Line(points={{20.5,91},
                {20.5,77.5},{7.74,77.5},{7.74,63.6}},             color={0,0,
                127}));
        connect(pump_type2_1.InFlow,tank. OutFlow) annotation (Line(points={{-90.86,
                -13.96},{-95.43,-13.96},{-95.43,-63},{-70.34,-63}},      color=
                {0,0,255}));
        connect(tank.InFlow,condenser. OutFlow_working) annotation (Line(points={{-46.4,
                -46.8},{-44,-46.8},{-44,-58.24},{-11.52,-58.24}},  color={0,0,
                255}));
        connect(turbine_simple.OutFlow,condenser. InFlow_working) annotation (
            Line(points={{54.6,-48.6},{54.6,-57.3},{36,-57.3},{36,-57.88}},
              color={0,0,255}));
        connect(hx1counter.OutFlow_working, turbine_simple.InFlow) annotation (
            Line(points={{11.58,9.56},{11.58,-27.22},{45,-27.22},{45,-37.8}},
              color={0,0,255}));
        connect(hx1counter.InFlow_working, pump_type2_1.OutFlow) annotation (
            Line(points={{-30,9.22},{-52,9.22},{-52,6.34},{-66.44,6.34}}, color=
               {0,0,255}));
        connect(condenser.OutFlow_liquid, outFlow_Tmflow.inFlow) annotation (
            Line(points={{36,-71.92},{59,-71.92},{59,-86.2},{70,-86.2}}, color=
                {28,108,200}));
        connect(inFlow_Tmflow.outFlow, condenser.InFlow_liquid) annotation (
            Line(points={{-72,-82},{-42,-82},{-42,-72.28},{-11.52,-72.28}},
              color={28,108,200}));
        connect(combiTimeTable.y[1], solarCollector1.I_in) annotation (Line(
              points={{-23.3,89},{-3.65,89},{-3.65,63.6},{-3.28,63.6}}, color={
                0,0,127}));
        connect(const3.y, valve.u) annotation (Line(points={{64.5,39},{64.25,39},
                {64.25,44},{76.64,44}}, color={0,0,127}));
        connect(mixer.outFlow, solarCollector1.inFlow) annotation (Line(points=
                {{-53.4,50},{-38,50},{-38,52},{-20,52}}, color={28,108,200}));
        connect(mixer.inFlow1, hx1counter.OutFlow_liquid) annotation (Line(
              points={{-59.4,44.8},{-69.7,44.8},{-69.7,22.48},{-30,22.48}},
              color={28,108,200}));
        connect(tank_liquid1.outFlow, mixer.inFlow2) annotation (Line(points={{
                27.2,70.4},{-59.4,70.4},{-59.4,55.8}}, color={28,108,200}));
        connect(tank_liquid1.inFlow, pump_liquid1.outFlow) annotation (Line(
              points={{38.2,80.2},{55.1,80.2},{55.1,72},{64,72}}, color={28,108,
                200}));
        connect(pump_liquid1.inFlow, valve.outFlow_2) annotation (Line(points={
                {84,72},{88,72},{88,44},{84.2,44}}, color={28,108,200}));
        connect(solarCollector1.outFlow, sens_liquid.inFlow) annotation (Line(
              points={{18,52.4},{30,52.4},{30,52.2},{40.4,52.2}}, color={28,108,
                200}));
        connect(sens_liquid.outFlow, valve.inFlow) annotation (Line(points={{
                51.2,52.2},{65.6,52.2},{65.6,50},{80,50}}, color={28,108,200}));
        connect(const4.y, add.u1) annotation (Line(points={{13.5,113},{43.75,
                113},{43.75,110},{49,110}}, color={0,0,127}));
        connect(sens_liquid.y, add.u2) annotation (Line(points={{45.6,55.2},{
                45.6,104.6},{49,104.6},{49,104}}, color={0,0,127}));
        connect(add.y, timeModel.u) annotation (Line(points={{60.5,107},{60.5,
                116.5},{72,116.5},{72,126.2}}, color={0,0,127}));
        connect(hx1counter.InFlow_liquid, pump_liquid.outFlow) annotation (Line(
              points={{11.58,22.82},{25.79,22.82},{25.79,14},{40,14}}, color={
                28,108,200}));
        connect(valve.outFlow_1, tank_liquid2.inFlow) annotation (Line(points={
                {79.88,38},{82,38},{82,0.2},{84.2,0.2}}, color={28,108,200}));
        connect(pump_liquid.inFlow, tank_liquid2.outFlow) annotation (Line(
              points={{60,14},{66,14},{66,-9.6},{73.2,-9.6}}, color={28,108,200}));
        connect(timeModel.y, gain.u) annotation (Line(points={{89,126},{96,126},
                {96,114},{117.2,114}}, color={0,0,127}));
        connect(const5.y, pump_liquid1.u) annotation (Line(points={{93.5,93},{
                83.75,93},{83.75,75.9},{73.9,75.9}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-78,-12},{8,-36}},
                lineColor={28,108,200},
                textString="需要保证冷凝器冷却足够，
使冷凝器出口工质达到过冷液态")}));
      end TestORCwithCollector;

      model TestNewORC_Colletor
        Components.main_components.Pump_liquid pump_liquid
          annotation (Placement(transformation(extent={{70,14},{50,34}})));
        Modelica.Blocks.Sources.Constant const(k=3.5)
          annotation (Placement(transformation(extent={{56,34},{46,44}})));
        Modelica.Blocks.Sources.Constant const1(k=800)
          annotation (Placement(transformation(extent={{-42,86},{-32,96}})));
        Modelica.Blocks.Sources.Constant const2(k=25)
          annotation (Placement(transformation(extent={{20,96},{30,106}})));
        Components.main_components.SolarCollector solarCollector1(
          K=1,
          N=20,
          Area=300,
          M=137,
          m_flow=4.5,
          Cp=2.34e3,
          T0=403.15)
          annotation (Placement(transformation(extent={{-10,42},{28,82}})));
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Uv=700,
          Cp_wall=0.39e3,
          max_drhode=200,
          witdth_x=0.1,
          Utp=5000,
          F_working=8.8,
          F_liquid=8.8,
          N=20,
          M_wall=10,
          M=300,
          Cp=2.34e3,
          m_liquid=4.5,
          Tin_start=312.15,
          Tout_start=388.15,
          T1_liquid=403.15,
          Tn_liquid=373.15,
          m_working=0.3)
          annotation (Placement(transformation(extent={{-20,8},{22,42}})));
        Components.main_components.Tank tank(redeclare package medium =
              ThermoCycle.Media.R245fa_CP, T0(displayUnit="degC") = 313.15)
          annotation (Placement(transformation(extent={{-70,-62},{-28,-26}})));
        Components.main_components.Pump_type2 pump_type2_1(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          v_s=1.2801e-5,
          T0=311.15,
          Np=35)
          annotation (Placement(transformation(extent={{-100,-26},{-34,32}})));
        Components.main_components.Hx1counter condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          N=10,
          Ul=1000,
          Uv=700,
          M=700,
          Cp=4.2e3,
          M_wall=50,
          Cp_wall=0.39e3,
          witdth_x=0.1,
          Utp=5000,
          F_working=13.8,
          F_liquid=13.8,
          max_drhode=200,
          m_liquid=5,
          m_working=0.3,
          Tin_start=371.15,
          Tout_start=303.15,
          T1_liquid=303.15,
          Tn_liquid=313.15,
          pstart=300000)
          annotation (Placement(transformation(extent={{46,-36},{-2,-72}})));
        Components.main_components.Turbine_simple turbine_simple(redeclare
            package medium = ThermoCycle.Media.R245fa_CP,
          n_turbine=0.7,
          V_s=3.121e-4,
          Np=38,
          Tin_start=383.15)
          annotation (Placement(transformation(extent={{50,-42},{70,-22}})));
        Components.source_sink_ports.OutFlow_Tmflow outFlow_Tmflow(m0=1, use_m0=
             false,
          N=10)
          annotation (Placement(transformation(extent={{80,-86},{100,-66}})));
        Components.source_sink_ports.InFlow_Tmflow inFlow_Tmflow(
          use_m0=true,
          T0=298.15,
          m0=-12)
          annotation (Placement(transformation(extent={{-82,-82},{-62,-62}})));
        Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(
          tableOnFile=true,
          columns={2},
          tableName="haha",
          fileName="haha.mat")
          annotation (Placement(transformation(extent={{-78,78},{-64,92}})));
        Components.main_components.sens_liquid sens_liquid
          annotation (Placement(transformation(extent={{46,52},{66,72}})));
        Components.main_components.Tank_liquid tank_liquid2(
          rou=825,
          Cp=2.34e3,
          Area=0.1,
          L0=0.1,
          T0=403.15)
          annotation (Placement(transformation(extent={{78,-4},{98,16}})));
        Components.Controls.SuperHeat superHeat(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{28,-14},{48,6}})));
      equation
        connect(const.y,pump_liquid. u) annotation (Line(points={{45.5,39},{
                44.5,39},{44.5,27.9},{59.9,27.9}},
                                               color={0,0,127}));
        connect(const2.y,solarCollector1. Tam_in) annotation (Line(points={{30.5,
                101},{30.5,87.5},{17.74,87.5},{17.74,73.6}},      color={0,0,
                127}));
        connect(pump_type2_1.InFlow,tank. OutFlow) annotation (Line(points={{-80.86,
                -3.96},{-85.43,-3.96},{-85.43,-53},{-60.34,-53}},        color=
                {0,0,255}));
        connect(tank.InFlow,condenser. OutFlow_working) annotation (Line(points={{-36.4,
                -36.8},{-34,-36.8},{-34,-48.24},{-1.52,-48.24}},   color={0,0,
                255}));
        connect(turbine_simple.OutFlow,condenser. InFlow_working) annotation (
            Line(points={{64.6,-38.6},{64.6,-47.3},{46,-47.3},{46,-47.88}},
              color={0,0,255}));
        connect(hx1counter.InFlow_working,pump_type2_1. OutFlow) annotation (
            Line(points={{-20,19.22},{-42,19.22},{-42,16.34},{-56.44,16.34}},
                                                                          color=
               {0,0,255}));
        connect(condenser.OutFlow_liquid, outFlow_Tmflow.inFlow) annotation (
            Line(points={{46,-61.92},{69,-61.92},{69,-76.2},{80,-76.2}}, color=
                {28,108,200}));
        connect(inFlow_Tmflow.outFlow, condenser.InFlow_liquid) annotation (
            Line(points={{-62,-72},{-32,-72},{-32,-62.28},{-1.52,-62.28}},
              color={28,108,200}));
        connect(solarCollector1.outFlow, sens_liquid.inFlow) annotation (Line(
              points={{28,62.4},{40,62.4},{40,62.2},{50.4,62.2}}, color={28,108,
                200}));
        connect(hx1counter.InFlow_liquid, pump_liquid.outFlow) annotation (Line(
              points={{21.58,32.82},{35.79,32.82},{35.79,24},{50,24}}, color={
                28,108,200}));
        connect(pump_liquid.inFlow, tank_liquid2.outFlow) annotation (Line(
              points={{70,24},{76,24},{76,0.4},{83.2,0.4}}, color={28,108,200}));
        connect(sens_liquid.outFlow, tank_liquid2.inFlow) annotation (Line(
              points={{61.2,62.2},{61.2,62.1},{94.2,62.1},{94.2,10.2}}, color={
                28,108,200}));
        connect(hx1counter.OutFlow_liquid, solarCollector1.inFlow) annotation (
            Line(points={{-20,32.48},{-20,62.24},{-10,62.24},{-10,62}}, color={
                28,108,200}));
        connect(hx1counter.OutFlow_working, superHeat.inFlow) annotation (Line(
              points={{21.58,19.56},{21.58,6.78},{28,6.78},{28,-3.8}}, color={0,
                0,255}));
        connect(superHeat.outFlow, turbine_simple.InFlow) annotation (Line(
              points={{48,-3.8},{52,-3.8},{52,-27.8},{55,-27.8}}, color={0,0,
                255}));
        connect(const1.y, solarCollector1.I_in) annotation (Line(points={{-31.5,
                91},{-31.5,89.5},{6.72,89.5},{6.72,73.6}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-64,-6},{22,-30}},
                lineColor={28,108,200},
                textString="晴天情况下系统循环
流量与响应关系")}));
      end TestNewORC_Colletor;

      model TestNewORC_Cloudy
        Components.main_components.Pump_liquid pump_liquid
          annotation (Placement(transformation(extent={{58,-6},{38,14}})));
        Modelica.Blocks.Sources.Constant const(k=3.5)
          annotation (Placement(transformation(extent={{44,14},{34,24}})));
        Modelica.Blocks.Sources.Constant const2(k=25)
          annotation (Placement(transformation(extent={{8,76},{18,86}})));
        Components.main_components.SolarCollector solarCollector1(
          K=1,
          N=20,
          Area=300,
          m_flow=4.5,
          Cp=2.34e3,
          M=137,
          T0=403.15)
          annotation (Placement(transformation(extent={{-22,22},{16,62}})));
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Uv=700,
          Cp_wall=0.39e3,
          witdth_x=0.1,
          Utp=5000,
          F_working=8.8,
          F_liquid=8.8,
          N=20,
          M_wall=10,
          M=300,
          Cp=2.34e3,
          m_liquid=4.5,
          m_working=0.3,
          max_drhode=100,
          V_working=1,
          V_liquid=1,
          Tin_start=312.15,
          Tout_start=388.15,
          T1_liquid=403.15,
          Tn_liquid=373.15)
          annotation (Placement(transformation(extent={{-32,-12},{10,22}})));
        Components.main_components.Tank tank(redeclare package medium =
              ThermoCycle.Media.R245fa_CP, T0(displayUnit="degC") = 313.15,
          Aera=1,
          L0=1)
          annotation (Placement(transformation(extent={{-82,-82},{-40,-46}})));
        Components.main_components.Pump_type2 pump_type2_1(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          v_s=1.2801e-5,
          T0=311.15,
          Np=40)
          annotation (Placement(transformation(extent={{-112,-46},{-46,12}})));
        Components.main_components.Hx1counter condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          N=10,
          Ul=1000,
          Uv=700,
          M=700,
          Cp=4.2e3,
          M_wall=50,
          Cp_wall=0.39e3,
          witdth_x=0.1,
          Utp=5000,
          F_working=13.8,
          F_liquid=13.8,
          m_liquid=5,
          m_working=0.3,
          max_drhode=100,
          Tin_start=371.15,
          Tout_start=303.15,
          T1_liquid=303.15,
          Tn_liquid=313.15,
          pstart=300000)
          annotation (Placement(transformation(extent={{34,-56},{-14,-92}})));
        Components.main_components.Turbine_simple turbine_simple(redeclare
            package medium = ThermoCycle.Media.R245fa_CP,
          n_turbine=0.7,
          V_s=3.121e-4,
          Np=40,
          Tin_start=383.15)
          annotation (Placement(transformation(extent={{38,-62},{58,-42}})));
        Components.source_sink_ports.OutFlow_Tmflow outFlow_Tmflow(m0=1, use_m0=
             false)
          annotation (Placement(transformation(extent={{68,-106},{88,-86}})));
        Components.source_sink_ports.InFlow_Tmflow inFlow_Tmflow(
          use_m0=true,
          T0=298.15,
          m0=-10)
          annotation (Placement(transformation(extent={{-94,-102},{-74,-82}})));
        Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(
          tableOnFile=true,
          columns={2},
          tableName="uaua",
          fileName="uaua.mat")
          annotation (Placement(transformation(extent={{-80,58},{-66,72}})));
        Components.main_components.sens_liquid sens_liquid
          annotation (Placement(transformation(extent={{34,32},{54,52}})));
        Components.main_components.Tank_liquid tank_liquid2(
          rou=825,
          Cp=2.34e3,
          Area=0.1,
          L0=0.1,
          T0=403.15)
          annotation (Placement(transformation(extent={{66,-24},{86,-4}})));
        Components.Controls.SuperHeat superHeat(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{16,-34},{36,-14}})));
        Modelica.Blocks.Sources.Constant const1(k=450)
          annotation (Placement(transformation(extent={{-64,82},{-54,92}})));
      equation
        connect(const.y,pump_liquid. u) annotation (Line(points={{33.5,19},{
                32.5,19},{32.5,7.9},{47.9,7.9}},
                                               color={0,0,127}));
        connect(const2.y,solarCollector1. Tam_in) annotation (Line(points={{18.5,81},
                {18.5,67.5},{5.74,67.5},{5.74,53.6}},             color={0,0,
                127}));
        connect(pump_type2_1.InFlow,tank. OutFlow) annotation (Line(points={{-92.86,
                -23.96},{-97.43,-23.96},{-97.43,-73},{-72.34,-73}},      color=
                {0,0,255}));
        connect(tank.InFlow,condenser. OutFlow_working) annotation (Line(points={{-48.4,
                -56.8},{-46,-56.8},{-46,-68.24},{-13.52,-68.24}},  color={0,0,
                255}));
        connect(turbine_simple.OutFlow,condenser. InFlow_working) annotation (
            Line(points={{52.6,-58.6},{52.6,-67.3},{34,-67.3},{34,-67.88}},
              color={0,0,255}));
        connect(hx1counter.InFlow_working,pump_type2_1. OutFlow) annotation (
            Line(points={{-32,-0.78},{-54,-0.78},{-54,-3.66},{-68.44,-3.66}},
                                                                          color=
               {0,0,255}));
        connect(condenser.OutFlow_liquid,outFlow_Tmflow. inFlow) annotation (
            Line(points={{34,-81.92},{57,-81.92},{57,-96.2},{68,-96.2}}, color=
                {28,108,200}));
        connect(inFlow_Tmflow.outFlow,condenser. InFlow_liquid) annotation (
            Line(points={{-74,-92},{-44,-92},{-44,-82.28},{-13.52,-82.28}},
              color={28,108,200}));
        connect(solarCollector1.outFlow,sens_liquid. inFlow) annotation (Line(
              points={{16,42.4},{28,42.4},{28,42.2},{38.4,42.2}}, color={28,108,
                200}));
        connect(hx1counter.InFlow_liquid,pump_liquid. outFlow) annotation (Line(
              points={{9.58,12.82},{23.79,12.82},{23.79,4},{38,4}},    color={
                28,108,200}));
        connect(pump_liquid.inFlow,tank_liquid2. outFlow) annotation (Line(
              points={{58,4},{64,4},{64,-19.6},{71.2,-19.6}},
                                                            color={28,108,200}));
        connect(sens_liquid.outFlow, tank_liquid2.inFlow) annotation (Line(
              points={{49.2,42.2},{49.2,42.1},{82.2,42.1},{82.2,-9.8}}, color={
                28,108,200}));
        connect(hx1counter.OutFlow_liquid, solarCollector1.inFlow) annotation (
            Line(points={{-32,12.48},{-32,42.24},{-22,42.24},{-22,42}}, color={
                28,108,200}));
        connect(hx1counter.OutFlow_working, superHeat.inFlow) annotation (Line(
              points={{9.58,-0.44},{9.58,-13.22},{16,-13.22},{16,-23.8}}, color=
               {0,0,255}));
        connect(superHeat.outFlow, turbine_simple.InFlow) annotation (Line(
              points={{36,-23.8},{40,-23.8},{40,-47.8},{43,-47.8}}, color={0,0,
                255}));
        connect(combiTimeTable.y[1], solarCollector1.I_in) annotation (Line(
              points={{-65.3,65},{-65.3,66.5},{-5.28,66.5},{-5.28,53.6}}, color=
               {0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-76,-26},{10,-50}},
                lineColor={28,108,200},
                textString="晴天情况下系统循环
流量与响应关系")}));
      end TestNewORC_Cloudy;

      model TestORCwithControl
        Components.main_components.Pump_liquid pump_liquid
          annotation (Placement(transformation(extent={{68,4},{48,24}})));
        Modelica.Blocks.Sources.Constant const(k=3.5)
          annotation (Placement(transformation(extent={{76,30},{66,40}})));
        Modelica.Blocks.Sources.Constant const2(k=25)
          annotation (Placement(transformation(extent={{18,86},{28,96}})));
        Components.main_components.SolarCollector solarCollector1(
          K=1,
          N=20,
          Area=300,
          m_flow=4.5,
          Cp=2.34e3,
          M=137,
          T0=403.15)
          annotation (Placement(transformation(extent={{-12,32},{26,72}})));
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Uv=700,
          Cp_wall=0.39e3,
          witdth_x=0.1,
          Utp=5000,
          F_working=8.8,
          F_liquid=8.8,
          N=20,
          M_wall=10,
          M=300,
          Cp=2.34e3,
          m_liquid=4.5,
          m_working=0.3,
          max_drhode=100,
          Tin_start=312.15,
          Tout_start=388.15,
          T1_liquid=403.15,
          Tn_liquid=373.15)
          annotation (Placement(transformation(extent={{-22,4},{20,38}})));
        Components.main_components.Tank tank(redeclare package medium =
              ThermoCycle.Media.R245fa_CP, T0(displayUnit="degC") = 313.15)
          annotation (Placement(transformation(extent={{-72,-72},{-30,-36}})));
        Components.main_components.Hx1counter condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          N=10,
          Ul=1000,
          Uv=700,
          M=700,
          Cp=4.2e3,
          M_wall=50,
          Cp_wall=0.39e3,
          witdth_x=0.1,
          Utp=5000,
          F_working=13.8,
          F_liquid=13.8,
          m_liquid=5,
          m_working=0.3,
          max_drhode=100,
          Tin_start=371.15,
          Tout_start=303.15,
          T1_liquid=303.15,
          Tn_liquid=313.15,
          pstart=300000)
          annotation (Placement(transformation(extent={{44,-46},{-4,-82}})));
        Components.source_sink_ports.OutFlow_Tmflow outFlow_Tmflow(m0=1, use_m0=
             false)
          annotation (Placement(transformation(extent={{78,-96},{98,-76}})));
        Components.source_sink_ports.InFlow_Tmflow inFlow_Tmflow(
          use_m0=true,
          m0=-10,
          T0=293.15)
          annotation (Placement(transformation(extent={{-84,-92},{-64,-72}})));
        Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(
          tableOnFile=true,
          columns={2},
          tableName="uaua",
          fileName="uaua.mat")
          annotation (Placement(transformation(extent={{-70,68},{-56,82}})));
        Components.main_components.sens_liquid sens_liquid
          annotation (Placement(transformation(extent={{44,42},{64,62}})));
        Components.main_components.Tank_liquid tank_liquid2(
          rou=825,
          Cp=2.34e3,
          Area=0.1,
          L0=0.1,
          T0=403.15)
          annotation (Placement(transformation(extent={{76,-14},{96,6}})));
        Components.Controls.SuperHeat superHeat(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{-6,-6},{6,6}},
              rotation=-90,
              origin={28,4})));
        Modelica.Blocks.Sources.Constant const1(k=450)
          annotation (Placement(transformation(extent={{-54,92},{-44,102}})));
        Components.main_components.pump_type3 pump_type3_1(
          v_s=1.2801e-5,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          T0=311.15)
          annotation (Placement(transformation(extent={{-42,-52},{-108,16}})));
        Modelica.Blocks.Sources.Constant const3(k=35)
          annotation (Placement(transformation(extent={{-118,16},{-108,26}})));
        Components.main_components.Turbine_simple_control
          turbine_simple_control(
          V_s=3.121e-4,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          n_turbine=0.7)
          annotation (Placement(transformation(extent={{48,-44},{68,-24}})));
        Components.Controls.sens_PT sens_PT(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{40,-10},{50,4}})));
        Components.Controls.control_unit1 control_unit1_1(
          P0=12,
          Y0=28,
          t=3000)
          annotation (Placement(transformation(extent={{32,-22},{46,-10}})));
        Modelica.Blocks.Math.Add add(k2=-1)
          annotation (Placement(transformation(extent={{2,-30},{14,-18}})));
        Modelica.Blocks.Sources.Constant const4(k=12)
          annotation (Placement(transformation(extent={{-20,-38},{-10,-28}})));
        Modelica.Blocks.Continuous.PI PI(initType=Modelica.Blocks.Types.Init.NoInit,
            T=1000)
          annotation (Placement(transformation(extent={{22,-38},{32,-28}})));
        Modelica.Blocks.Math.Add add1(k2=+1)
          annotation (Placement(transformation(extent={{112,-28},{124,-16}})));
        Modelica.Blocks.Sources.Constant const5(k=40)
          annotation (Placement(transformation(extent={{68,-60},{78,-50}})));
      equation
        connect(const.y,pump_liquid. u) annotation (Line(points={{65.5,35},{
                42.5,35},{42.5,17.9},{57.9,17.9}},
                                               color={0,0,127}));
        connect(const2.y,solarCollector1. Tam_in) annotation (Line(points={{28.5,91},
                {28.5,77.5},{15.74,77.5},{15.74,63.6}},           color={0,0,
                127}));
        connect(tank.InFlow,condenser. OutFlow_working) annotation (Line(points={{-38.4,
                -46.8},{-36,-46.8},{-36,-58.24},{-3.52,-58.24}},   color={0,0,
                255}));
        connect(condenser.OutFlow_liquid,outFlow_Tmflow. inFlow) annotation (
            Line(points={{44,-71.92},{67,-71.92},{67,-86.2},{78,-86.2}}, color=
                {28,108,200}));
        connect(inFlow_Tmflow.outFlow,condenser. InFlow_liquid) annotation (
            Line(points={{-64,-82},{-34,-82},{-34,-72.28},{-3.52,-72.28}},
              color={28,108,200}));
        connect(solarCollector1.outFlow,sens_liquid. inFlow) annotation (Line(
              points={{26,52.4},{38,52.4},{38,52.2},{48.4,52.2}}, color={28,108,
                200}));
        connect(hx1counter.InFlow_liquid,pump_liquid. outFlow) annotation (Line(
              points={{19.58,28.82},{33.79,28.82},{33.79,14},{48,14}}, color={
                28,108,200}));
        connect(pump_liquid.inFlow,tank_liquid2. outFlow) annotation (Line(
              points={{68,14},{74,14},{74,-9.6},{81.2,-9.6}},
                                                            color={28,108,200}));
        connect(sens_liquid.outFlow,tank_liquid2. inFlow) annotation (Line(
              points={{59.2,52.2},{59.2,52.1},{92.2,52.1},{92.2,0.2}},  color={
                28,108,200}));
        connect(hx1counter.OutFlow_liquid,solarCollector1. inFlow) annotation (
            Line(points={{-22,28.48},{-22,52.24},{-12,52.24},{-12,52}}, color={
                28,108,200}));
        connect(hx1counter.OutFlow_working,superHeat. inFlow) annotation (Line(
              points={{19.58,15.56},{19.58,16.78},{28.12,16.78},{28.12,10}},
                                                                          color=
               {0,0,255}));
        connect(combiTimeTable.y[1],solarCollector1. I_in) annotation (Line(
              points={{-55.3,75},{-55.3,76.5},{4.72,76.5},{4.72,63.6}},   color=
               {0,0,127}));
        connect(pump_type3_1.InFlow, tank.OutFlow) annotation (Line(points={{-61.14,
                -26.16},{-61.14,-37.08},{-62.34,-37.08},{-62.34,-63}},
              color={0,0,255}));
        connect(pump_type3_1.OutFlow, hx1counter.InFlow_working) annotation (
            Line(points={{-85.56,-2.36},{-51.78,-2.36},{-51.78,15.22},{-22,
                15.22}},
              color={0,0,255}));
        connect(condenser.InFlow_working, turbine_simple_control.OutFlow)
          annotation (Line(points={{44,-57.88},{53,-57.88},{53,-40.6},{62.6,
                -40.6}}, color={0,0,255}));
        connect(const3.y, pump_type3_1.u) annotation (Line(points={{-107.5,21},
                {-96.75,21},{-96.75,-13.92},{-84.9,-13.92}}, color={0,0,127}));
        connect(superHeat.outFlow, sens_PT.inFlow) annotation (Line(points={{
                28.12,-2},{30,-2},{30,-2.86},{40,-2.86}}, color={0,0,255}));
        connect(sens_PT.outFlow, turbine_simple_control.InFlow) annotation (
            Line(points={{50,-2.86},{53,-2.86},{53,-29.8}}, color={0,0,255}));
        connect(sens_PT.y_P, control_unit1_1.u) annotation (Line(points={{42.6,
                0.92},{42.6,-7.54},{33.26,-7.54},{33.26,-16}}, color={0,0,127}));
        connect(sens_PT.y_P, add.u1) annotation (Line(points={{42.6,0.92},{
                -20.7,0.92},{-20.7,-20.4},{0.8,-20.4}}, color={0,0,127}));
        connect(const4.y, add.u2) annotation (Line(points={{-9.5,-33},{-4.75,
                -33},{-4.75,-27.6},{0.8,-27.6}}, color={0,0,127}));
        connect(add.y, PI.u) annotation (Line(points={{14.6,-24},{20,-24},{20,
                -33},{21,-33}}, color={0,0,127}));
        connect(control_unit1_1.y, add1.u1) annotation (Line(points={{43.62,-16},
                {70,-16},{70,-18.4},{110.8,-18.4}},color={0,0,127}));
        connect(PI.y, add1.u2) annotation (Line(points={{32.5,-33},{47.25,-33},
                {47.25,-25.6},{110.8,-25.6}},color={0,0,127}));
        connect(const5.y, turbine_simple_control.u) annotation (Line(points={{
                78.5,-55},{78.5,-34.5},{62.3,-34.5},{62.3,-33.3}}, color={0,0,
                127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestORCwithControl;

      model TestORCwithControl2
        Components.main_components.Pump_liquid pump_liquid
          annotation (Placement(transformation(extent={{78,14},{58,34}})));
        Modelica.Blocks.Sources.Constant const(k=3.5)
          annotation (Placement(transformation(extent={{86,40},{76,50}})));
        Modelica.Blocks.Sources.Constant const2(k=25)
          annotation (Placement(transformation(extent={{28,96},{38,106}})));
        Components.main_components.SolarCollector solarCollector1(
          K=1,
          N=20,
          Area=300,
          m_flow=4.5,
          Cp=2.34e3,
          M=137,
          T0=403.15)
          annotation (Placement(transformation(extent={{-2,42},{36,82}})));
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Uv=700,
          Cp_wall=0.39e3,
          witdth_x=0.1,
          Utp=5000,
          F_working=8.8,
          F_liquid=8.8,
          N=20,
          M_wall=10,
          M=300,
          Cp=2.34e3,
          m_liquid=4.5,
          m_working=0.3,
          max_drhode=100,
          Tin_start=312.15,
          Tout_start=388.15,
          T1_liquid=403.15,
          Tn_liquid=373.15)
          annotation (Placement(transformation(extent={{-12,14},{30,48}})));
        Components.main_components.Tank tank(redeclare package medium =
              ThermoCycle.Media.R245fa_CP, T0(displayUnit="degC") = 313.15)
          annotation (Placement(transformation(extent={{-62,-62},{-20,-26}})));
        Components.main_components.Hx1counter condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Uv=700,
          M=700,
          Cp=4.2e3,
          M_wall=50,
          Cp_wall=0.39e3,
          witdth_x=0.1,
          Utp=5000,
          F_working=13.8,
          F_liquid=13.8,
          m_liquid=5,
          m_working=0.3,
          max_drhode=100,
          N=20,
          Tin_start=371.15,
          Tout_start=303.15,
          T1_liquid=303.15,
          Tn_liquid=313.15,
          pstart=300000)
          annotation (Placement(transformation(extent={{54,-36},{6,-72}})));
        Components.source_sink_ports.OutFlow_Tmflow outFlow_Tmflow(m0=1, use_m0=
             false,
          N=20)
          annotation (Placement(transformation(extent={{88,-86},{108,-66}})));
        Components.source_sink_ports.InFlow_Tmflow inFlow_Tmflow(
          use_m0=true,
          m0=-18,
          T0=298.15)
          annotation (Placement(transformation(extent={{-74,-82},{-54,-62}})));
        Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(
          tableOnFile=true,
          columns={2},
          tableName="uaua",
          fileName="uaua.mat")
          annotation (Placement(transformation(extent={{-60,78},{-46,92}})));
        Components.main_components.sens_liquid sens_liquid
          annotation (Placement(transformation(extent={{54,52},{74,72}})));
        Components.main_components.Tank_liquid tank_liquid2(
          rou=825,
          Cp=2.34e3,
          Area=0.1,
          L0=0.1,
          T0=403.15)
          annotation (Placement(transformation(extent={{98,8},{118,28}})));
        Components.Controls.SuperHeat superHeat(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{-6,-6},{6,6}},
              rotation=-90,
              origin={38,18})));
        Modelica.Blocks.Sources.Constant const1(k=30)
          annotation (Placement(transformation(extent={{-178,-14},{-170,-6}})));
        Components.main_components.pump_type3 pump_type3_1(
          v_s=1.2801e-5,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          T0=311.15)
          annotation (Placement(transformation(extent={{-28,-36},{-94,32}})));
        Modelica.Blocks.Sources.Constant const3(k=28)
          annotation (Placement(transformation(extent={{-106,36},{-96,46}})));
        Components.Controls.controlsource controlsource
          annotation (Placement(transformation(extent={{-96,48},{-76,68}})));
        Components.main_components.Turbine_simple_control
          turbine_simple_control(
          V_s=3.121e-4,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          n_turbine=0.7)
          annotation (Placement(transformation(extent={{64,-40},{84,-20}})));
        Components.Controls.sens_PT sens_PT(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{50,0},{60,14}})));
        Components.Controls.control_unit1 control_unit1_1(
          P0=12,
          Y0=28,
          t=3000)
          annotation (Placement(transformation(extent={{38,-14},{52,-2}})));
        Modelica.Blocks.Math.Add add(k2=-1)
          annotation (Placement(transformation(extent={{12,-20},{24,-8}})));
        Modelica.Blocks.Sources.Constant const4(k=12)
          annotation (Placement(transformation(extent={{-10,-28},{0,-18}})));
        Modelica.Blocks.Continuous.PI PI(initType=Modelica.Blocks.Types.Init.NoInit,
            T=1000)
          annotation (Placement(transformation(extent={{36,-32},{46,-22}})));
        Modelica.Blocks.Math.Add add1(k2=+1)
          annotation (Placement(transformation(extent={{94,-20},{106,-8}})));
        Modelica.Blocks.Math.Add add2(k2=-1)
          annotation (Placement(transformation(extent={{-154,-10},{-142,2}})));
        Modelica.Blocks.Math.Gain gain(k=1/7)
          annotation (Placement(transformation(extent={{-134,0},{-128,6}})));
        Modelica.Blocks.Math.Add add3(k2=+1)
          annotation (Placement(transformation(extent={{-112,0},{-100,12}})));
        Modelica.Blocks.Continuous.PI PI1(
          initType=Modelica.Blocks.Types.Init.NoInit,
          T=8000,
          k=0.65) annotation (Placement(transformation(extent={{-118,-30},{-108,
                  -20}})));
        Modelica.Blocks.Math.Add add4(k2=+1)
          annotation (Placement(transformation(extent={{-98,-28},{-86,-16}})));
      equation
        connect(const.y,pump_liquid. u) annotation (Line(points={{75.5,45},{
                52.5,45},{52.5,27.9},{67.9,27.9}},
                                               color={0,0,127}));
        connect(const2.y,solarCollector1. Tam_in) annotation (Line(points={{38.5,
                101},{38.5,87.5},{25.74,87.5},{25.74,73.6}},      color={0,0,
                127}));
        connect(tank.InFlow,condenser. OutFlow_working) annotation (Line(points={{-28.4,
                -36.8},{-26,-36.8},{-26,-48.24},{6.48,-48.24}},    color={0,0,
                255}));
        connect(condenser.OutFlow_liquid,outFlow_Tmflow. inFlow) annotation (
            Line(points={{54,-61.92},{77,-61.92},{77,-76.2},{88,-76.2}}, color=
                {28,108,200}));
        connect(inFlow_Tmflow.outFlow,condenser. InFlow_liquid) annotation (
            Line(points={{-54,-72},{-24,-72},{-24,-62.28},{6.48,-62.28}},
              color={28,108,200}));
        connect(solarCollector1.outFlow,sens_liquid. inFlow) annotation (Line(
              points={{36,62.4},{48,62.4},{48,62.2},{58.4,62.2}}, color={28,108,
                200}));
        connect(hx1counter.InFlow_liquid,pump_liquid. outFlow) annotation (Line(
              points={{29.58,38.82},{43.79,38.82},{43.79,24},{58,24}}, color={
                28,108,200}));
        connect(pump_liquid.inFlow,tank_liquid2. outFlow) annotation (Line(
              points={{78,24},{84,24},{84,12.4},{103.2,12.4}},
                                                            color={28,108,200}));
        connect(sens_liquid.outFlow,tank_liquid2. inFlow) annotation (Line(
              points={{69.2,62.2},{69.2,62.1},{114.2,62.1},{114.2,22.2}},
                                                                        color={
                28,108,200}));
        connect(hx1counter.OutFlow_liquid,solarCollector1. inFlow) annotation (
            Line(points={{-12,38.48},{-12,62.24},{-2,62.24},{-2,62}},   color={
                28,108,200}));
        connect(hx1counter.OutFlow_working,superHeat. inFlow) annotation (Line(
              points={{29.58,25.56},{29.58,26.78},{38.12,26.78},{38.12,24}},
                                                                          color=
               {0,0,255}));
        connect(combiTimeTable.y[1],solarCollector1. I_in) annotation (Line(
              points={{-45.3,85},{-45.3,86.5},{14.72,86.5},{14.72,73.6}}, color=
               {0,0,127}));
        connect(pump_type3_1.InFlow,tank. OutFlow) annotation (Line(points={{-47.14,
                -10.16},{-47.14,-27.08},{-52.34,-27.08},{-52.34,-53}},
              color={0,0,255}));
        connect(pump_type3_1.OutFlow,hx1counter. InFlow_working) annotation (
            Line(points={{-71.56,13.64},{-41.78,13.64},{-41.78,25.22},{-12,
                25.22}},
              color={0,0,255}));
        connect(condenser.InFlow_working,turbine_simple_control. OutFlow)
          annotation (Line(points={{54,-47.88},{79,-47.88},{79,-36.6},{78.6,
                -36.6}}, color={0,0,255}));
        connect(superHeat.outFlow, sens_PT.inFlow) annotation (Line(points={{
                38.12,12},{40,12},{40,7.14},{50,7.14}}, color={0,0,255}));
        connect(sens_PT.outFlow, turbine_simple_control.InFlow) annotation (
            Line(points={{60,7.14},{69,7.14},{69,-25.8}}, color={0,0,255}));
        connect(sens_PT.y_P, control_unit1_1.u) annotation (Line(points={{52.6,
                10.92},{52.6,4.46},{39.26,4.46},{39.26,-8}}, color={0,0,127}));
        connect(sens_PT.y_P, add.u1) annotation (Line(points={{52.6,10.92},{
                -4.7,10.92},{-4.7,-10.4},{10.8,-10.4}}, color={0,0,127}));
        connect(const4.y, add.u2) annotation (Line(points={{0.5,-23},{5.25,-23},
                {5.25,-17.6},{10.8,-17.6}}, color={0,0,127}));
        connect(add.y, PI.u) annotation (Line(points={{24.6,-14},{30,-14},{30,
                -27},{35,-27}}, color={0,0,127}));
        connect(control_unit1_1.y, add1.u1) annotation (Line(points={{49.62,-8},
                {80,-8},{80,-10.4},{92.8,-10.4}}, color={0,0,127}));
        connect(PI.y, add1.u2) annotation (Line(points={{46.5,-27},{53.25,-27},
                {53.25,-17.6},{92.8,-17.6}}, color={0,0,127}));
        connect(add1.y, turbine_simple_control.u) annotation (Line(points={{106.6,
                -14},{126,-14},{126,-29.3},{78.3,-29.3}},       color={0,0,127}));
        connect(superHeat.y, add2.u1) annotation (Line(points={{39.92,17.88},{
                -162.04,17.88},{-162.04,-0.4},{-155.2,-0.4}}, color={0,0,127}));
        connect(const1.y, add2.u2) annotation (Line(points={{-169.6,-10},{
                -166.75,-10},{-166.75,-7.6},{-155.2,-7.6}}, color={0,0,127}));
        connect(add2.y, gain.u) annotation (Line(points={{-141.4,-4},{-138,-4},
                {-138,3},{-134.6,3}}, color={0,0,127}));
        connect(const3.y, add3.u1) annotation (Line(points={{-95.5,41},{-95.5,
                20.5},{-113.2,20.5},{-113.2,9.6}},
                                                 color={0,0,127}));
        connect(gain.y, add3.u2) annotation (Line(points={{-127.7,3},{-114,3},{
                -114,2.4},{-113.2,2.4}}, color={0,0,127}));
        connect(add2.y, PI1.u) annotation (Line(points={{-141.4,-4},{-130,-4},{
                -130,-25},{-119,-25}}, color={0,0,127}));
        connect(add3.y, add4.u1) annotation (Line(points={{-99.4,6},{-90,6},{
                -90,-18.4},{-99.2,-18.4}}, color={0,0,127}));
        connect(PI1.y, add4.u2) annotation (Line(points={{-107.5,-25},{-96.75,
                -25},{-96.75,-25.6},{-99.2,-25.6}}, color={0,0,127}));
        connect(add4.y, pump_type3_1.u) annotation (Line(points={{-85.4,-22},{
                -80,-22},{-80,2.08},{-70.9,2.08}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestORCwithControl2;

      model TestORCwithControl2_Submodel
        Components.main_components.Pump_liquid pump_liquid
          annotation (Placement(transformation(extent={{88,24},{68,44}})));
        Modelica.Blocks.Sources.Constant const(k=3.5)
          annotation (Placement(transformation(extent={{96,50},{86,60}})));
        Modelica.Blocks.Sources.Constant const2(k=25)
          annotation (Placement(transformation(extent={{74,92},{64,102}})));
        Components.main_components.SolarCollector solarCollector1(
          K=1,
          N=20,
          Area=300,
          m_flow=4.5,
          Cp=2.34e3,
          M=137,
          T0=403.15)
          annotation (Placement(transformation(extent={{8,52},{46,92}})));
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Uv=700,
          Cp_wall=0.39e3,
          witdth_x=0.1,
          Utp=5000,
          F_working=8.8,
          F_liquid=8.8,
          N=20,
          M_wall=10,
          M=300,
          Cp=2.34e3,
          m_liquid=4.5,
          m_working=0.3,
          max_drhode=100,
          Tin_start=312.15,
          Tout_start=388.15,
          T1_liquid=403.15,
          Tn_liquid=373.15)
          annotation (Placement(transformation(extent={{-2,24},{40,58}})));
        Components.main_components.Tank tank(redeclare package medium =
              ThermoCycle.Media.R245fa_CP, T0(displayUnit="degC") = 313.15)
          annotation (Placement(transformation(extent={{-52,-52},{-10,-16}})));
        Components.main_components.Hx1counter condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Uv=700,
          M=700,
          Cp=4.2e3,
          M_wall=50,
          Cp_wall=0.39e3,
          witdth_x=0.1,
          Utp=5000,
          F_working=13.8,
          F_liquid=13.8,
          m_liquid=5,
          m_working=0.3,
          max_drhode=100,
          N=20,
          Tin_start=371.15,
          Tout_start=303.15,
          T1_liquid=303.15,
          Tn_liquid=313.15,
          pstart=300000)
          annotation (Placement(transformation(extent={{64,-26},{16,-62}})));
        Components.source_sink_ports.OutFlow_Tmflow outFlow_Tmflow(m0=1, use_m0=
             false)
          annotation (Placement(transformation(extent={{98,-76},{118,-56}})));
        Components.source_sink_ports.InFlow_Tmflow inFlow_Tmflow(
          use_m0=true,
          m0=-18,
          T0=298.15)
          annotation (Placement(transformation(extent={{-64,-72},{-44,-52}})));
        Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(
          tableOnFile=true,
          columns={2},
          tableName="uaua",
          fileName="uaua.mat")
          annotation (Placement(transformation(extent={{-50,88},{-36,102}})));
        Components.main_components.sens_liquid sens_liquid
          annotation (Placement(transformation(extent={{64,62},{84,82}})));
        Components.main_components.Tank_liquid tank_liquid2(
          rou=825,
          Cp=2.34e3,
          Area=0.1,
          L0=0.1,
          T0=403.15)
          annotation (Placement(transformation(extent={{98,30},{118,50}})));
        Components.Controls.SuperHeat superHeat(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{-6,-6},{6,6}},
              rotation=-90,
              origin={48,28})));
        Components.main_components.pump_type3 pump_type3_1(
          v_s=1.2801e-5,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          T0=311.15)
          annotation (Placement(transformation(extent={{-18,-26},{-84,42}})));
        Modelica.Blocks.Sources.Constant const3(k=28)
          annotation (Placement(transformation(extent={{-96,-30},{-86,-20}})));
        Components.main_components.Turbine_simple_control
          turbine_simple_control(
          V_s=3.121e-4,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          n_turbine=0.7)
          annotation (Placement(transformation(extent={{60,-28},{84,4}})));
        Components.Controls.sens_PT sens_PT(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{60,10},{70,24}})));
        Modelica.Blocks.Math.Add add(k2=-1)
          annotation (Placement(transformation(extent={{102,2},{114,-10}})));
        Modelica.Blocks.Sources.Constant const4(k=12)
          annotation (Placement(transformation(extent={{120,10},{110,20}})));
        control1 control1_1 annotation (Placement(transformation(rotation=0,
                extent={{128,-10},{148,10}})));
        control2 control2_1 annotation (Placement(transformation(rotation=0,
                extent={{-92,22},{-72,2}})));
      equation
        connect(const.y,pump_liquid. u) annotation (Line(points={{85.5,55},{
                62.5,55},{62.5,37.9},{77.9,37.9}},
                                               color={0,0,127}));
        connect(const2.y,solarCollector1. Tam_in) annotation (Line(points={{63.5,97},
                {63.5,97.5},{35.74,97.5},{35.74,83.6}},           color={0,0,
                127}));
        connect(tank.InFlow,condenser. OutFlow_working) annotation (Line(points={{-18.4,
                -26.8},{-16,-26.8},{-16,-38.24},{16.48,-38.24}},   color={0,0,
                255}));
        connect(condenser.OutFlow_liquid,outFlow_Tmflow. inFlow) annotation (
            Line(points={{64,-51.92},{87,-51.92},{87,-66.2},{98,-66.2}}, color=
                {28,108,200}));
        connect(inFlow_Tmflow.outFlow,condenser. InFlow_liquid) annotation (
            Line(points={{-44,-62},{-14,-62},{-14,-52.28},{16.48,-52.28}},
              color={28,108,200}));
        connect(solarCollector1.outFlow,sens_liquid. inFlow) annotation (Line(
              points={{46,72.4},{58,72.4},{58,72.2},{68.4,72.2}}, color={28,108,
                200}));
        connect(hx1counter.InFlow_liquid,pump_liquid. outFlow) annotation (Line(
              points={{39.58,48.82},{53.79,48.82},{53.79,34},{68,34}}, color={
                28,108,200}));
        connect(pump_liquid.inFlow,tank_liquid2. outFlow) annotation (Line(
              points={{88,34},{94,34},{94,34.4},{103.2,34.4}},
                                                            color={28,108,200}));
        connect(sens_liquid.outFlow,tank_liquid2. inFlow) annotation (Line(
              points={{79.2,72.2},{79.2,72.1},{114.2,72.1},{114.2,44.2}},
                                                                        color={
                28,108,200}));
        connect(hx1counter.OutFlow_liquid,solarCollector1. inFlow) annotation (
            Line(points={{-2,48.48},{-2,72.24},{8,72.24},{8,72}},       color={
                28,108,200}));
        connect(hx1counter.OutFlow_working,superHeat. inFlow) annotation (Line(
              points={{39.58,35.56},{39.58,36.78},{48.12,36.78},{48.12,34}},
                                                                          color=
               {0,0,255}));
        connect(combiTimeTable.y[1],solarCollector1. I_in) annotation (Line(
              points={{-35.3,95},{-35.3,96.5},{24.72,96.5},{24.72,83.6}}, color=
               {0,0,127}));
        connect(pump_type3_1.InFlow,tank. OutFlow) annotation (Line(points={{-37.14,
                -0.16},{-37.14,-17.08},{-42.34,-17.08},{-42.34,-43}},
              color={0,0,255}));
        connect(pump_type3_1.OutFlow,hx1counter. InFlow_working) annotation (
            Line(points={{-61.56,23.64},{-31.78,23.64},{-31.78,35.22},{-2,35.22}},
              color={0,0,255}));
        connect(condenser.InFlow_working,turbine_simple_control. OutFlow)
          annotation (Line(points={{64,-37.88},{73,-37.88},{73,-22.56},{77.52,
                -22.56}},color={0,0,255}));
        connect(superHeat.outFlow,sens_PT. inFlow) annotation (Line(points={{48.12,
                22},{50,22},{50,17.14},{60,17.14}},     color={0,0,255}));
        connect(sens_PT.outFlow,turbine_simple_control. InFlow) annotation (
            Line(points={{70,17.14},{70,-5.28},{66,-5.28}},
                                                          color={0,0,255}));
        connect(sens_PT.y_P,add. u1) annotation (Line(points={{62.6,20.92},{
                93.3,20.92},{93.3,-7.6},{100.8,-7.6}},  color={0,0,127}));
        connect(const4.y,add. u2) annotation (Line(points={{109.5,15},{99.25,15},
                {99.25,-0.4},{100.8,-0.4}}, color={0,0,127}));
        connect(add.y, control1_1.u) annotation (Line(points={{114.6,-4},{122,
                -4},{128,-4}}, color={0,0,127}));
        connect(superHeat.y, control2_1.u1) annotation (Line(points={{49.92,
                27.88},{-92,27.88},{-92,16.2}}, color={0,0,127}));
        connect(control2_1.y, pump_type3_1.u) annotation (Line(points={{-71.8,
                11.8},{-70,11.8},{-70,12.08},{-60.9,12.08}}, color={0,0,127}));
        connect(control1_1.y, turbine_simple_control.u) annotation (Line(points=
               {{148.2,1.2},{156,1.2},{156,-20},{86,-20},{86,-10.88},{77.16,
                -10.88}}, color={0,0,127}));
        connect(const3.y, control2_1.u2) annotation (Line(points={{-85.5,-25},{
                -66,-25},{-66,-4},{-96,-4},{-96,7.8},{-91.8,7.8}}, color={0,0,
                127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestORCwithControl2_Submodel;

      model control1
        Components.Controls.control_unit1 control_unit1_1(
          P0=12,
          Y0=28,
          t=3000)
          annotation (Placement(transformation(extent={{142,-4},{156,8}})));
        Modelica.Blocks.Continuous.PI PI(initType=Modelica.Blocks.Types.Init.NoInit,
            T=1000)
          annotation (Placement(transformation(extent={{124,-38},{134,-28}})));
        Modelica.Blocks.Math.Add add1(k2=+1)
          annotation (Placement(transformation(extent={{156,-36},{168,-24}})));
        Modelica.Blocks.Interfaces.RealInput u annotation (Placement(
              transformation(rotation=0, extent={{-110,-50},{-90,-30}}),
              iconTransformation(extent={{-110,-50},{-90,-30}})));
        Modelica.Blocks.Interfaces.RealOutput y annotation (Placement(
              transformation(rotation=0, extent={{92,2},{112,22}}),
              iconTransformation(extent={{92,2},{112,22}})));
      equation
        connect(control_unit1_1.y,add1. u1) annotation (Line(points={{153.62,2},
                {188,2},{188,-26.4},{154.8,-26.4}},
                                                  color={0,0,127}));
        connect(PI.y,add1. u2) annotation (Line(points={{134.5,-33},{133.25,-33},
                {133.25,-33.6},{154.8,-33.6}},
                                             color={0,0,127}));
        connect(u, PI.u) annotation (Line(points={{-100,-40},{18,-40},{18,-34},
                {70,-34},{70,-33},{123,-33}}, color={0,0,127}));
        connect(y, add1.y) annotation (Line(points={{102,12},{136,12},{136,-30},
                {168.6,-30}}, color={0,0,127}));
        connect(u, control_unit1_1.u) annotation (Line(points={{-100,-40},{22,
                -40},{22,2},{143.26,2}}, color={0,0,127}));
        annotation (Icon(graphics={Rectangle(
                extent={{-100,74},{100,-78}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid), Text(
                extent={{-76,40},{60,-32}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid,
                textString="control1")}));
      end control1;

      model control2
        Modelica.Blocks.Sources.Constant const1(k=30)
          annotation (Placement(transformation(extent={{-168,-4},{-160,4}})));
        Modelica.Blocks.Math.Add add2(k2=-1)
          annotation (Placement(transformation(extent={{-144,0},{-132,12}})));
        Modelica.Blocks.Math.Gain gain(k=1/7)
          annotation (Placement(transformation(extent={{-124,12},{-118,18}})));
        Modelica.Blocks.Math.Add add3(k2=+1)
          annotation (Placement(transformation(extent={{-106,6},{-94,18}})));
        Modelica.Blocks.Continuous.PI PI1(
          initType=Modelica.Blocks.Types.Init.NoInit,
          T=8000,
          k=0.65)
          annotation (Placement(transformation(extent={{-108,-20},{-98,-10}})));
        Modelica.Blocks.Math.Add add4(k2=+1)
          annotation (Placement(transformation(extent={{-88,-18},{-76,-6}})));
        Modelica.Blocks.Interfaces.RealInput u1 annotation (Placement(
              transformation(rotation=0, extent={{-110,-52},{-90,-32}}),
              iconTransformation(extent={{-110,-52},{-90,-32}})));
        Modelica.Blocks.Interfaces.RealInput u2 annotation (Placement(
              transformation(rotation=0, extent={{-108,32},{-88,52}}),
              iconTransformation(extent={{-108,32},{-88,52}})));
        Modelica.Blocks.Interfaces.RealOutput y annotation (Placement(
              transformation(rotation=0, extent={{92,-8},{112,12}}),
              iconTransformation(extent={{92,-8},{112,12}})));
      equation
        connect(const1.y,add2. u2) annotation (Line(points={{-159.6,0},{-156.75,
                0},{-156.75,2.4},{-145.2,2.4}},             color={0,0,127}));
        connect(add2.y,gain. u) annotation (Line(points={{-131.4,6},{-128,6},{
                -128,15},{-124.6,15}},color={0,0,127}));
        connect(gain.y,add3. u2) annotation (Line(points={{-117.7,15},{-104,15},
                {-104,8.4},{-107.2,8.4}},color={0,0,127}));
        connect(add2.y,PI1. u) annotation (Line(points={{-131.4,6},{-120,6},{
                -120,-15},{-109,-15}}, color={0,0,127}));
        connect(add3.y,add4. u1) annotation (Line(points={{-93.4,12},{-76,12},{
                -76,-8.4},{-89.2,-8.4}},   color={0,0,127}));
        connect(PI1.y,add4. u2) annotation (Line(points={{-97.5,-15},{-86.75,
                -15},{-86.75,-15.6},{-89.2,-15.6}}, color={0,0,127}));
        connect(u1, add2.u1) annotation (Line(points={{-100,-42},{-100,9.6},{
                -145.2,9.6}}, color={0,0,127}));
        connect(u2, add3.u1) annotation (Line(points={{-98,42},{-107.2,42},{
                -107.2,58},{-107.2,15.6}}, color={0,0,127}));
        connect(y, add4.y) annotation (Line(points={{102,2},{12,2},{12,-12},{
                -75.4,-12}}, color={0,0,127}));
        annotation (Icon(graphics={Rectangle(
                extent={{-100,72},{100,-68}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid), Text(
                extent={{-80,38},{62,-34}},
                lineColor={28,108,200},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid,
                textString="Control2")}));
      end control2;

      model TestPTC_ORC_sunny_day_stepchange
        Components.main_components.Pump_liquid pump_liquid
          annotation (Placement(transformation(extent={{64,4},{44,24}})));
        Modelica.Blocks.Sources.Constant const(k=3.5)
          annotation (Placement(transformation(extent={{50,24},{40,34}})));
        Modelica.Blocks.Sources.Constant const1(k=800)
          annotation (Placement(transformation(extent={{-48,76},{-38,86}})));
        Modelica.Blocks.Sources.Constant const2(k=25)
          annotation (Placement(transformation(extent={{14,86},{24,96}})));
        Components.main_components.SolarCollector solarCollector1(
          K=1,
          N=20,
          Area=300,
          M=137,
          m_flow=4.5,
          Cp=2.34e3,
          T0=403.15)
          annotation (Placement(transformation(extent={{-16,32},{22,72}})));
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Uv=700,
          Cp_wall=0.39e3,
          max_drhode=200,
          witdth_x=0.1,
          Utp=5000,
          F_working=8.8,
          F_liquid=8.8,
          N=20,
          M_wall=10,
          M=300,
          Cp=2.34e3,
          m_liquid=4.5,
          Tin_start=312.15,
          Tout_start=388.15,
          T1_liquid=403.15,
          Tn_liquid=373.15,
          m_working=0.3)
          annotation (Placement(transformation(extent={{-26,-2},{16,32}})));
        Components.main_components.Tank tank(redeclare package medium =
              ThermoCycle.Media.R245fa_CP, T0(displayUnit="degC") = 313.15,
          L0=0.5)
          annotation (Placement(transformation(extent={{-76,-72},{-34,-36}})));
        Components.main_components.Hx1counter condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          N=10,
          Ul=1000,
          Uv=700,
          M=700,
          Cp=4.2e3,
          M_wall=50,
          Cp_wall=0.39e3,
          witdth_x=0.1,
          Utp=5000,
          F_working=13.8,
          F_liquid=13.8,
          max_drhode=200,
          m_liquid=5,
          m_working=0.3,
          Tin_start=371.15,
          Tout_start=303.15,
          T1_liquid=303.15,
          Tn_liquid=313.15,
          pstart=300000)
          annotation (Placement(transformation(extent={{40,-46},{-8,-82}})));
        Components.main_components.Turbine_simple turbine_simple(redeclare
            package medium = ThermoCycle.Media.R245fa_CP,
          n_turbine=0.7,
          V_s=3.121e-4,
          Np=38,
          Tin_start=383.15)
          annotation (Placement(transformation(extent={{44,-52},{64,-32}})));
        Components.source_sink_ports.OutFlow_Tmflow outFlow_Tmflow(m0=1, use_m0=
             false,
          N=10)
          annotation (Placement(transformation(extent={{74,-96},{94,-76}})));
        Components.source_sink_ports.InFlow_Tmflow inFlow_Tmflow(
          use_m0=true,
          T0=298.15,
          m0=-15)
          annotation (Placement(transformation(extent={{-88,-92},{-68,-72}})));
        Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(
          tableOnFile=true,
          columns={2},
          tableName="haha",
          fileName="haha.mat")
          annotation (Placement(transformation(extent={{-80,68},{-66,82}})));
        Components.main_components.sens_liquid sens_liquid
          annotation (Placement(transformation(extent={{40,42},{60,62}})));
        Components.main_components.Tank_liquid tank_liquid2(
          rou=825,
          Cp=2.34e3,
          Area=0.1,
          L0=0.1,
          T0=403.15)
          annotation (Placement(transformation(extent={{72,-14},{92,6}})));
        Components.Controls.SuperHeat superHeat(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{22,-24},{42,-4}})));
        Components.main_components.pump_type3 pump_type3_1(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          T0=311.15,
          v_s=1.2801e-5)
          annotation (Placement(transformation(extent={{-58,-58},{-100,-6}})));
        Modelica.Blocks.Sources.Constant const3(k=35)
          annotation (Placement(transformation(extent={{-124,-4},{-114,6}})));
        Modelica.Blocks.Sources.Step step(
          startTime=1000,
          height=10,
          offset=35) annotation (Placement(transformation(extent={{-118,24},{
                  -102,40}})));
      equation
        connect(const.y,pump_liquid. u) annotation (Line(points={{39.5,29},{
                38.5,29},{38.5,17.9},{53.9,17.9}},
                                               color={0,0,127}));
        connect(const2.y,solarCollector1. Tam_in) annotation (Line(points={{24.5,91},
                {24.5,77.5},{11.74,77.5},{11.74,63.6}},           color={0,0,
                127}));
        connect(tank.InFlow,condenser. OutFlow_working) annotation (Line(points={{-42.4,
                -46.8},{-40,-46.8},{-40,-58.24},{-7.52,-58.24}},   color={0,0,
                255}));
        connect(turbine_simple.OutFlow,condenser. InFlow_working) annotation (
            Line(points={{58.6,-48.6},{58.6,-57.3},{40,-57.3},{40,-57.88}},
              color={0,0,255}));
        connect(condenser.OutFlow_liquid,outFlow_Tmflow. inFlow) annotation (
            Line(points={{40,-71.92},{63,-71.92},{63,-86.2},{74,-86.2}}, color=
                {28,108,200}));
        connect(inFlow_Tmflow.outFlow,condenser. InFlow_liquid) annotation (
            Line(points={{-68,-82},{-38,-82},{-38,-72.28},{-7.52,-72.28}},
              color={28,108,200}));
        connect(solarCollector1.outFlow,sens_liquid. inFlow) annotation (Line(
              points={{22,52.4},{34,52.4},{34,52.2},{44.4,52.2}}, color={28,108,
                200}));
        connect(hx1counter.InFlow_liquid,pump_liquid. outFlow) annotation (Line(
              points={{15.58,22.82},{29.79,22.82},{29.79,14},{44,14}}, color={
                28,108,200}));
        connect(pump_liquid.inFlow,tank_liquid2. outFlow) annotation (Line(
              points={{64,14},{70,14},{70,-9.6},{77.2,-9.6}},
                                                            color={28,108,200}));
        connect(sens_liquid.outFlow,tank_liquid2. inFlow) annotation (Line(
              points={{55.2,52.2},{55.2,52.1},{88.2,52.1},{88.2,0.2}},  color={
                28,108,200}));
        connect(hx1counter.OutFlow_liquid,solarCollector1. inFlow) annotation (
            Line(points={{-26,22.48},{-26,52.24},{-16,52.24},{-16,52}}, color={
                28,108,200}));
        connect(hx1counter.OutFlow_working,superHeat. inFlow) annotation (Line(
              points={{15.58,9.56},{15.58,-3.22},{22,-3.22},{22,-13.8}},
                                                                       color={0,
                0,255}));
        connect(superHeat.outFlow,turbine_simple. InFlow) annotation (Line(
              points={{42,-13.8},{46,-13.8},{46,-37.8},{49,-37.8}},
                                                                  color={0,0,
                255}));
        connect(const1.y, solarCollector1.I_in) annotation (Line(points={{-37.5,
                81},{-37.5,79.5},{0.72,79.5},{0.72,63.6}}, color={0,0,127}));
        connect(pump_type3_1.OutFlow, hx1counter.InFlow_working) annotation (
            Line(points={{-85.72,-20.04},{-90.14,-20.04},{-90.14,9.22},{-26,
                9.22}}, color={0,0,255}));
        connect(pump_type3_1.InFlow, tank.OutFlow) annotation (Line(points={{
                -70.18,-38.24},{-70.18,-62.12},{-66.34,-62.12},{-66.34,-63}},
              color={0,0,255}));
        connect(step.y, pump_type3_1.u) annotation (Line(points={{-101.2,32},{
                -94,32},{-94,-28.88},{-85.3,-28.88}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-70,-16},{16,-40}},
                lineColor={28,108,200},
                textString="晴天情况下系统循环
流量与响应关系")}));
      end TestPTC_ORC_sunny_day_stepchange;

      model TestPTC_ORC_sunny_day_stepchange2
        Modelica.Blocks.Sources.Step step(
          startTime=1000,
          offset=38,
          height=-10)
                     annotation (Placement(transformation(extent={{58,-8},{74,8}})));
        Components.main_components.Turbine_simple_control
          turbine_simple_control(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          n_turbine=0.7,
          V_s=3.121e-4,
          Tin_start=383.15)
          annotation (Placement(transformation(extent={{14,-28},{34,-8}})));
        Components.main_components.Pump_liquid pump_liquid
          annotation (Placement(transformation(extent={{24,20},{4,40}})));
        Modelica.Blocks.Sources.Constant const(k=3.5)
          annotation (Placement(transformation(extent={{10,40},{0,50}})));
        Modelica.Blocks.Sources.Constant const1(k=800)
          annotation (Placement(transformation(extent={{-84,100},{-74,110}})));
        Modelica.Blocks.Sources.Constant const2(k=25)
          annotation (Placement(transformation(extent={{-26,102},{-16,112}})));
        Components.main_components.SolarCollector solarCollector1(
          K=1,
          N=20,
          Area=300,
          M=137,
          m_flow=4.5,
          Cp=2.34e3,
          T0=403.15)
          annotation (Placement(transformation(extent={{-56,48},{-18,88}})));
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Uv=700,
          Cp_wall=0.39e3,
          max_drhode=200,
          witdth_x=0.1,
          Utp=5000,
          F_working=8.8,
          F_liquid=8.8,
          N=20,
          M_wall=10,
          M=300,
          Cp=2.34e3,
          m_liquid=4.5,
          Tin_start=312.15,
          Tout_start=388.15,
          T1_liquid=403.15,
          Tn_liquid=373.15,
          m_working=0.3)
          annotation (Placement(transformation(extent={{-66,14},{-24,48}})));
        Components.main_components.Tank tank(redeclare package medium =
              ThermoCycle.Media.R245fa_CP, T0(displayUnit="degC") = 313.15,
          L0=0.5)
          annotation (Placement(transformation(extent={{-116,-56},{-74,-20}})));
        Components.main_components.Hx1counter condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          N=10,
          Ul=1000,
          Uv=700,
          M=700,
          Cp=4.2e3,
          M_wall=50,
          Cp_wall=0.39e3,
          witdth_x=0.1,
          Utp=5000,
          F_working=13.8,
          F_liquid=13.8,
          max_drhode=200,
          m_liquid=5,
          m_working=0.3,
          Tin_start=371.15,
          Tout_start=303.15,
          T1_liquid=303.15,
          Tn_liquid=313.15,
          pstart=300000)
          annotation (Placement(transformation(extent={{0,-30},{-48,-66}})));
        Components.source_sink_ports.OutFlow_Tmflow outFlow_Tmflow(m0=1, use_m0=
             false)
          annotation (Placement(transformation(extent={{34,-80},{54,-60}})));
        Components.source_sink_ports.InFlow_Tmflow inFlow_Tmflow(
          use_m0=true,
          T0=298.15,
          m0=-15)
          annotation (Placement(transformation(extent={{-128,-76},{-108,-56}})));
        Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(
          tableOnFile=true,
          columns={2},
          tableName="haha",
          fileName="haha.mat")
          annotation (Placement(transformation(extent={{-164,100},{-150,114}})));
        Components.main_components.sens_liquid sens_liquid
          annotation (Placement(transformation(extent={{0,58},{20,78}})));
        Components.main_components.Tank_liquid tank_liquid2(
          rou=825,
          Cp=2.34e3,
          Area=0.1,
          L0=0.5,
          T0=403.15)
          annotation (Placement(transformation(extent={{32,2},{52,22}})));
        Components.Controls.SuperHeat superHeat(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{-18,-8},{2,12}})));
        Components.main_components.pump_type3 pump_type3_1(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          T0=311.15,
          v_s=1.2801e-5)
          annotation (Placement(transformation(extent={{-98,-42},{-140,10}})));
        Modelica.Blocks.Sources.Constant const3(k=45)
          annotation (Placement(transformation(extent={{-164,12},{-154,22}})));
        Modelica.Blocks.Sources.Step step1(
          offset=0,
          height=-670,
          startTime=10000)
                     annotation (Placement(transformation(extent={{-160,82},{
                  -152,90}})));
        Modelica.Blocks.Sources.Constant const4(k=45)
          annotation (Placement(transformation(extent={{32,-46},{42,-36}})));
        Modelica.Blocks.Sources.Step step2(height=670, startTime=10300)
          annotation (Placement(transformation(extent={{-162,54},{-152,64}})));
        Modelica.Blocks.Math.Add3 add3_1
          annotation (Placement(transformation(extent={{-116,80},{-96,100}})));
      equation
        connect(const.y,pump_liquid. u) annotation (Line(points={{-0.5,45},{
                -1.5,45},{-1.5,33.9},{13.9,33.9}},
                                               color={0,0,127}));
        connect(const2.y,solarCollector1. Tam_in) annotation (Line(points={{-15.5,
                107},{-15.5,93.5},{-28.26,93.5},{-28.26,79.6}},   color={0,0,
                127}));
        connect(tank.InFlow,condenser. OutFlow_working) annotation (Line(points={{-82.4,
                -30.8},{-80,-30.8},{-80,-42.24},{-47.52,-42.24}},  color={0,0,
                255}));
        connect(condenser.OutFlow_liquid,outFlow_Tmflow. inFlow) annotation (
            Line(points={{0,-55.92},{23,-55.92},{23,-70.2},{34,-70.2}},  color=
                {28,108,200}));
        connect(inFlow_Tmflow.outFlow,condenser. InFlow_liquid) annotation (
            Line(points={{-108,-66},{-78,-66},{-78,-56.28},{-47.52,-56.28}},
              color={28,108,200}));
        connect(solarCollector1.outFlow,sens_liquid. inFlow) annotation (Line(
              points={{-18,68.4},{-6,68.4},{-6,68.2},{4.4,68.2}}, color={28,108,
                200}));
        connect(hx1counter.InFlow_liquid,pump_liquid. outFlow) annotation (Line(
              points={{-24.42,38.82},{-10.21,38.82},{-10.21,30},{4,30}},
                                                                       color={
                28,108,200}));
        connect(pump_liquid.inFlow,tank_liquid2. outFlow) annotation (Line(
              points={{24,30},{30,30},{30,6.4},{37.2,6.4}}, color={28,108,200}));
        connect(sens_liquid.outFlow,tank_liquid2. inFlow) annotation (Line(
              points={{15.2,68.2},{15.2,68.1},{48.2,68.1},{48.2,16.2}}, color={
                28,108,200}));
        connect(hx1counter.OutFlow_liquid,solarCollector1. inFlow) annotation (
            Line(points={{-66,38.48},{-66,68.24},{-56,68.24},{-56,68}}, color={
                28,108,200}));
        connect(hx1counter.OutFlow_working,superHeat. inFlow) annotation (Line(
              points={{-24.42,25.56},{-24.42,12.78},{-18,12.78},{-18,2.2}},
                                                                       color={0,
                0,255}));
        connect(pump_type3_1.OutFlow,hx1counter. InFlow_working) annotation (
            Line(points={{-125.72,-4.04},{-130.14,-4.04},{-130.14,25.22},{-66,
                25.22}},color={0,0,255}));
        connect(pump_type3_1.InFlow,tank. OutFlow) annotation (Line(points={{-110.18,
                -22.24},{-110.18,-46.12},{-106.34,-46.12},{-106.34,-47}},
              color={0,0,255}));
        connect(const3.y, pump_type3_1.u) annotation (Line(points={{-153.5,17},
                {-153.5,-12.5},{-125.3,-12.5},{-125.3,-12.88}}, color={0,0,127}));
        connect(superHeat.outFlow, turbine_simple_control.InFlow) annotation (
            Line(points={{2,2.2},{12,2.2},{12,-13.8},{19,-13.8}}, color={0,0,
                255}));
        connect(condenser.InFlow_working, turbine_simple_control.OutFlow)
          annotation (Line(points={{0,-41.88},{16,-41.88},{16,-24.6},{28.6,
                -24.6}}, color={0,0,255}));
        connect(const4.y, turbine_simple_control.u) annotation (Line(points={{
                42.5,-41},{42.5,-17.5},{28.3,-17.5},{28.3,-17.3}}, color={0,0,
                127}));
        connect(combiTimeTable.y[1], add3_1.u1) annotation (Line(points={{
                -149.3,107},{-135.65,107},{-135.65,98},{-118,98}}, color={0,0,
                127}));
        connect(step1.y, add3_1.u2) annotation (Line(points={{-151.6,86},{-136,
                86},{-136,90},{-118,90}}, color={0,0,127}));
        connect(step2.y, add3_1.u3) annotation (Line(points={{-151.5,59},{
                -136.75,59},{-136.75,82},{-118,82}}, color={0,0,127}));
        connect(add3_1.y, solarCollector1.I_in) annotation (Line(points={{-95,
                90},{-39.28,90},{-39.28,79.6}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-110,0},{-24,-24}},
                lineColor={28,108,200},
                textString="晴天情况下系统循环
流量与响应关系")}));
      end TestPTC_ORC_sunny_day_stepchange2;

      model TestORCwithControl_pump
        Modelica.Blocks.Math.Add add(k2=-1)
          annotation (Placement(transformation(extent={{-18,-20},{-6,-8}})));
        Modelica.Blocks.Sources.Constant const4(k=20)
          annotation (Placement(transformation(extent={{-52,-12},{-42,-2}})));
        Modelica.Blocks.Continuous.PID PID(
          initType=Modelica.Blocks.Types.InitPID.NoInit,
          xi_start=0,
          xd_start=0,
          Td=100,
          k=0.02,
          Ti=150)
          annotation (Placement(transformation(extent={{-4,-38},{-16,-26}})));
        Components.main_components.Pump_liquid pump_liquid
          annotation (Placement(transformation(extent={{68,4},{48,24}})));
        Modelica.Blocks.Sources.Constant const(k=3.5)
          annotation (Placement(transformation(extent={{54,24},{44,34}})));
        Modelica.Blocks.Sources.Constant const2(k=25)
          annotation (Placement(transformation(extent={{18,86},{28,96}})));
        Components.main_components.SolarCollector solarCollector1(
          K=1,
          N=20,
          Area=300,
          m_flow=4.5,
          Cp=2.34e3,
          M=137,
          T0=403.15)
          annotation (Placement(transformation(extent={{-12,32},{26,72}})));
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Uv=700,
          Cp_wall=0.39e3,
          witdth_x=0.1,
          Utp=5000,
          F_working=8.8,
          F_liquid=8.8,
          N=20,
          M_wall=10,
          M=300,
          Cp=2.34e3,
          m_liquid=4.5,
          m_working=0.3,
          max_drhode=100,
          V_working=1,
          V_liquid=1,
          Tin_start=312.15,
          Tout_start=388.15,
          T1_liquid=403.15,
          Tn_liquid=373.15)
          annotation (Placement(transformation(extent={{-22,-2},{20,32}})));
        Components.main_components.Tank tank(redeclare package medium =
              ThermoCycle.Media.R245fa_CP, T0(displayUnit="degC") = 313.15,
          Aera=1,
          L0=1)
          annotation (Placement(transformation(extent={{-72,-72},{-30,-36}})));
        Components.main_components.Hx1counter condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          N=10,
          Ul=1000,
          Uv=700,
          M=700,
          Cp=4.2e3,
          M_wall=50,
          Cp_wall=0.39e3,
          witdth_x=0.1,
          Utp=5000,
          F_working=13.8,
          F_liquid=13.8,
          m_liquid=5,
          m_working=0.3,
          max_drhode=100,
          Tin_start=371.15,
          Tout_start=303.15,
          T1_liquid=303.15,
          Tn_liquid=313.15,
          pstart=300000)
          annotation (Placement(transformation(extent={{44,-46},{-4,-82}})));
        Components.main_components.Turbine_simple turbine_simple(redeclare
            package medium = ThermoCycle.Media.R245fa_CP,
          n_turbine=0.7,
          V_s=3.121e-4,
          Np=40,
          Tin_start=383.15)
          annotation (Placement(transformation(extent={{48,-52},{68,-32}})));
        Components.source_sink_ports.OutFlow_Tmflow outFlow_Tmflow(m0=1, use_m0=
             false)
          annotation (Placement(transformation(extent={{78,-96},{98,-76}})));
        Components.source_sink_ports.InFlow_Tmflow inFlow_Tmflow(
          use_m0=true,
          m0=-10,
          T0=293.15)
          annotation (Placement(transformation(extent={{-84,-92},{-64,-72}})));
        Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(
          tableOnFile=true,
          columns={2},
          tableName="uaua",
          fileName="uaua.mat")
          annotation (Placement(transformation(extent={{-70,68},{-56,82}})));
        Components.main_components.sens_liquid sens_liquid
          annotation (Placement(transformation(extent={{44,42},{64,62}})));
        Components.main_components.Tank_liquid tank_liquid2(
          rou=825,
          Cp=2.34e3,
          Area=0.1,
          L0=0.1,
          T0=403.15)
          annotation (Placement(transformation(extent={{76,-14},{96,6}})));
        Components.Controls.SuperHeat superHeat(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{26,-14},{46,6}})));
        Modelica.Blocks.Sources.Constant const1(k=450)
          annotation (Placement(transformation(extent={{-54,92},{-44,102}})));
        Components.main_components.pump_type3 pump_type3_1(
          v_s=1.2801e-5,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          T0=311.15)
          annotation (Placement(transformation(extent={{-124,-46},{-58,22}})));
        Modelica.Blocks.Math.Add add1(k2=+1)
          annotation (Placement(transformation(extent={{-38,-36},{-50,-24}})));
        Modelica.Blocks.Sources.Constant const3(k=40)
          annotation (Placement(transformation(extent={{-34,-48},{-24,-38}})));
        Modelica.Blocks.Sources.Constant const5(k=35)
          annotation (Placement(transformation(extent={{-92,-42},{-82,-32}})));
      equation
        connect(const4.y,add. u2) annotation (Line(points={{-41.5,-7},{-30.75,
                -7},{-30.75,-17.6},{-19.2,-17.6}},
                                                 color={0,0,127}));
        connect(add.y, PID.u) annotation (Line(points={{-5.4,-14},{0,-14},{0,
                -32},{-2.8,-32}}, color={0,0,127}));
        connect(const.y,pump_liquid. u) annotation (Line(points={{43.5,29},{
                42.5,29},{42.5,17.9},{57.9,17.9}},
                                               color={0,0,127}));
        connect(const2.y,solarCollector1. Tam_in) annotation (Line(points={{28.5,91},
                {28.5,77.5},{15.74,77.5},{15.74,63.6}},           color={0,0,
                127}));
        connect(tank.InFlow,condenser. OutFlow_working) annotation (Line(points={{-38.4,
                -46.8},{-36,-46.8},{-36,-58.24},{-3.52,-58.24}},   color={0,0,
                255}));
        connect(turbine_simple.OutFlow,condenser. InFlow_working) annotation (
            Line(points={{62.6,-48.6},{62.6,-57.3},{44,-57.3},{44,-57.88}},
              color={0,0,255}));
        connect(condenser.OutFlow_liquid,outFlow_Tmflow. inFlow) annotation (
            Line(points={{44,-71.92},{67,-71.92},{67,-86.2},{78,-86.2}}, color=
                {28,108,200}));
        connect(inFlow_Tmflow.outFlow,condenser. InFlow_liquid) annotation (
            Line(points={{-64,-82},{-34,-82},{-34,-72.28},{-3.52,-72.28}},
              color={28,108,200}));
        connect(solarCollector1.outFlow,sens_liquid. inFlow) annotation (Line(
              points={{26,52.4},{38,52.4},{38,52.2},{48.4,52.2}}, color={28,108,
                200}));
        connect(hx1counter.InFlow_liquid,pump_liquid. outFlow) annotation (Line(
              points={{19.58,22.82},{33.79,22.82},{33.79,14},{48,14}}, color={
                28,108,200}));
        connect(pump_liquid.inFlow,tank_liquid2. outFlow) annotation (Line(
              points={{68,14},{74,14},{74,-9.6},{81.2,-9.6}},
                                                            color={28,108,200}));
        connect(sens_liquid.outFlow,tank_liquid2. inFlow) annotation (Line(
              points={{59.2,52.2},{59.2,52.1},{92.2,52.1},{92.2,0.2}},  color={
                28,108,200}));
        connect(hx1counter.OutFlow_liquid,solarCollector1. inFlow) annotation (
            Line(points={{-22,22.48},{-22,52.24},{-12,52.24},{-12,52}}, color={
                28,108,200}));
        connect(hx1counter.OutFlow_working,superHeat. inFlow) annotation (Line(
              points={{19.58,9.56},{19.58,-3.22},{26,-3.22},{26,-3.8}},   color=
               {0,0,255}));
        connect(superHeat.outFlow,turbine_simple. InFlow) annotation (Line(
              points={{46,-3.8},{50,-3.8},{50,-37.8},{53,-37.8}},   color={0,0,
                255}));
        connect(combiTimeTable.y[1],solarCollector1. I_in) annotation (Line(
              points={{-55.3,75},{-55.3,76.5},{4.72,76.5},{4.72,63.6}},   color=
               {0,0,127}));
        connect(superHeat.y, add.u1) annotation (Line(points={{36.2,-0.8},{
                -23.9,-0.8},{-23.9,-10.4},{-19.2,-10.4}}, color={0,0,127}));
        connect(pump_type3_1.InFlow, tank.OutFlow) annotation (Line(points={{
                -104.86,-20.16},{-104.86,-64.08},{-62.34,-64.08},{-62.34,-63}},
              color={0,0,255}));
        connect(pump_type3_1.OutFlow, hx1counter.InFlow_working) annotation (
            Line(points={{-80.44,3.64},{-51.22,3.64},{-51.22,9.22},{-22,9.22}},
              color={0,0,255}));
        connect(PID.y, add1.u2) annotation (Line(points={{-16.6,-32},{-36.8,-32},
                {-36.8,-33.6}}, color={0,0,127}));
        connect(const3.y, add1.u1) annotation (Line(points={{-23.5,-43},{-21.75,
                -43},{-21.75,-26.4},{-36.8,-26.4}}, color={0,0,127}));
        connect(pump_type3_1.u, add1.y) annotation (Line(points={{-81.1,-7.92},
                {-65.55,-7.92},{-65.55,-30},{-50.6,-30}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestORCwithControl_pump;
    end TestCycle;

    package TestCycle2
      model CycleOne_ICE
        Components.main_components.Hx1counter hx1counter
          annotation (Placement(transformation(extent={{36,-50},{-4,-80}})));
        Components.main_components.Hx1counter hx1counter1 annotation (Placement(
              transformation(
              extent={{17,20},{-17,-20}},
              rotation=-90,
              origin={-50,-23})));
        Components.main_components.Hx1counter hx1counter2
          annotation (Placement(transformation(extent={{-28,2},{4,30}})));
        Components.main_components.Turbine_simple turbine_simple
          annotation (Placement(transformation(extent={{30,-30},{60,4}})));
        Components.main_components.Hx1counter_liquid hx1counter_liquid(
          N=10,
          Cp1=2.34e3,
          U1=1300,
          U2=800,
          F_liquid1=10,
          F_liquid2=10,
          Cp2=1.15e3,
          Cp_wall=400)
          annotation (Placement(transformation(extent={{-26,42},{10,68}})));
        Components.main_components.Pump_liquid pump_liquid(m0=0.1)
          annotation (Placement(transformation(extent={{48,16},{28,36}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid(T0(displayUnit=
               "K") = 808.15)
          annotation (Placement(transformation(extent={{52,64},{32,84}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid(m0=0.2981)
          annotation (Placement(transformation(extent={{-52,66},{-72,86}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1
          annotation (Placement(transformation(extent={{48,-84},{68,-64}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid1
          annotation (Placement(transformation(extent={{-72,-86},{-52,-66}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid2
          annotation (Placement(transformation(extent={{-92,-4},{-72,16}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid2
          annotation (Placement(transformation(extent={{-76,-46},{-96,-26}})));
        Components.main_components.Tank tank
          annotation (Placement(transformation(extent={{-32,-66},{-12,-46}})));
        Components.main_components.Pump_simple pump_simple
          annotation (Placement(transformation(extent={{-60,-66},{-32,-42}})));
      equation
        connect(hx1counter1.OutFlow_working, hx1counter2.InFlow_working)
          annotation (Line(points={{-43.6,-6.34},{-43.6,-4.17},{-28,-4.17},{-28,
                11.24}}, color={0,0,255}));
        connect(hx1counter2.OutFlow_working, turbine_simple.InFlow) annotation (
           Line(points={{3.68,11.52},{20.84,11.52},{20.84,-5.86},{37.5,-5.86}},
              color={0,0,255}));
        connect(pump_liquid.outFlow, hx1counter2.InFlow_liquid) annotation (
            Line(points={{28,26},{16,26},{16,22.44},{3.68,22.44}}, color={28,
                108,200}));
        connect(hx1counter_liquid.outflow1, pump_liquid.inFlow) annotation (
            Line(points={{9.64,49.8},{56.82,49.8},{56.82,26},{48,26}}, color={
                28,108,200}));
        connect(hx1counter2.OutFlow_liquid, hx1counter_liquid.inflow1)
          annotation (Line(points={{-28,22.16},{-46,22.16},{-46,49.8},{-26.36,
                49.8}}, color={28,108,200}));
        connect(inFlow_liquid.A, hx1counter_liquid.inflow2) annotation (Line(
              points={{32.6,73.8},{30.7,73.8},{30.7,61.24},{9.64,61.24}}, color=
               {28,108,200}));
        connect(outFlow_liquid.A, hx1counter_liquid.outflow2) annotation (Line(
              points={{-52,76},{-40,76},{-40,61.5},{-26,61.5}}, color={28,108,
                200}));
        connect(hx1counter.OutFlow_liquid, outFlow_liquid1.A) annotation (Line(
              points={{36,-71.6},{38,-71.6},{38,-74},{48,-74}}, color={28,108,
                200}));
        connect(inFlow_liquid1.A, hx1counter.InFlow_liquid) annotation (Line(
              points={{-52.6,-76.2},{-30.3,-76.2},{-30.3,-71.9},{-3.6,-71.9}},
              color={28,108,200}));
        connect(inFlow_liquid2.A, hx1counter1.InFlow_liquid) annotation (Line(
              points={{-72.6,5.8},{-72.6,4.9},{-59.2,4.9},{-59.2,-6.34}}, color=
               {28,108,200}));
        connect(outFlow_liquid2.A, hx1counter1.OutFlow_liquid) annotation (Line(
              points={{-76,-36},{-70,-36},{-70,-40},{-58.8,-40}}, color={28,108,
                200}));
        connect(turbine_simple.OutFlow, hx1counter.InFlow_working) annotation (
            Line(points={{51.9,-24.22},{51.9,-42.11},{36,-42.11},{36,-59.9}},
              color={0,0,255}));
        connect(tank.InFlow, hx1counter.OutFlow_working) annotation (Line(
              points={{-16,-52},{-12,-52},{-12,-60.2},{-3.6,-60.2}}, color={0,0,
                255}));
        connect(pump_simple.InFlow, tank.OutFlow) annotation (Line(points={{
                -51.88,-56.88},{-51.94,-56.88},{-51.94,-61},{-27.4,-61}}, color=
               {0,0,255}));
        connect(pump_simple.OutFlow, hx1counter1.InFlow_working) annotation (
            Line(points={{-41.52,-48.48},{-41.52,-45.24},{-43.2,-45.24},{-43.2,
                -40}}, color={0,0,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-28,-8},{32,-48}},
                lineColor={28,108,200},
                fillColor={170,213,255},
                fillPattern=FillPattern.Forward,
                textString="内燃机余热回收系统")}));
      end CycleOne_ICE;

      model CycleOneICE_part1
        Components.main_components.Hx1counter_liquid Gas_oil_Hx(
          Cp1=2.34e3,
          U1=1300,
          U2=800,
          Cp_wall=400,
          M_wall=20,
          F_liquid1=3,
          F_liquid2=3,
          V_liquid1=0.7,
          V_liquid2=0.7,
          M1=500,
          N=20,
          T1_liquid(displayUnit="degC") = 447.15,
          T1n_liquid(displayUnit="degC") = 499.55,
          T2_liquid(displayUnit="degC") = 808.15,
          T2n_liquid(displayUnit="degC") = 455.42,
          Cp2=1.109e3,
          M2=0.5)
          annotation (Placement(transformation(extent={{-24,66},{12,92}})));
        Components.main_components.Hx1counter Evaporater(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Utp=3000,
          Uv=800,
          U=1000,
          M=500,
          Cp=2.34e3,
          M_wall=20,
          Cp_wall=400,
          F_working=2.3,
          F_liquid=2.3,
          m_liquid=1,
          N=20,
          m_working=0.7,
          Tin_start=365.33,
          Tout_start=393.15,
          T1_liquid=499.55,
          Tn_liquid=447.85,
          pstart=1895000)
          annotation (Placement(transformation(extent={{-14,14},{18,42}})));
        Components.main_components.Pump_liquid pump_liquid
          annotation (Placement(transformation(extent={{50,26},{30,46}})));
        Modelica.Blocks.Sources.Constant const(k=1)
          annotation (Placement(transformation(extent={{16,46},{26,56}})));
        Components.Controls.sens_PT sens_PT(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{32,8},{52,28}})));
        Components.Controls.sens_PT sens_PT1(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{-46,14},{-26,34}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(m0=2.1)
          annotation (Placement(transformation(extent={{-74,-26},{-94,-6}})));
        Components.main_components.Hx1counter preHeater(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Utp=3000,
          Uv=800,
          U=1000,
          M_wall=20,
          Cp_wall=400,
          N=20,
          M=700,
          Cp=4.2e3,
          m_liquid=2.1,
          F_working=10,
          F_liquid=10,
          m_working=0.7,
          Tin_start=296.45,
          Tout_start=365.35,
          T1_liquid=366.15,
          Tn_liquid=358.45,
          pstart=1895000) annotation (Placement(transformation(
              extent={{-10,-11},{10,11}},
              rotation=90,
              origin={-63,8})));
        Components.main_components.Turbine_simple turbine_simple(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          Np=35,
          V_s=1.904e-4,
          pin_start=1895000,
          pout_start=176000,
          Tin_start=398.15)
          annotation (Placement(transformation(extent={{76,-30},{104,4}})));
        Components.main_components.Hx1counter Condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Utp=3000,
          Uv=800,
          U=1000,
          M_wall=20,
          Cp_wall=400,
          N=20,
          M=700,
          Cp=4.2e3,
          m_working=0.7,
          m_liquid=10,
          F_working=20,
          F_liquid=20,
          Tin_start=347.11,
          Tout_start=295.38,
          T1_liquid=290.15,
          Tn_liquid=294.15,
          pstart=176000) annotation (Placement(transformation(
              extent={{-16.5,-14.5},{16.5,14.5}},
              rotation=180,
              origin={23.5,-46.5})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid2(m0=10)
          annotation (Placement(transformation(extent={{66,-74},{86,-54}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid2(T0(
              displayUnit="degC") = 290.15)
          annotation (Placement(transformation(extent={{32,-90},{12,-70}})));
        Components.Controls.sens_PT sens_PT2(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{2,-52},{-18,-32}})));
        Components.main_components.Pump_type2 pump_type2_1(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          v_s=1.0395e-5,
          Np=50,
          n_pump=0.7,
          p0=176000,
          pstart_out=1895000,
          T0=295.15)
          annotation (Placement(transformation(extent={{-56,-58},{-82,-30}})));
        Components.main_components.Tank tank(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          Aera=0.065,
          L0=0.65,
          p0=176000,
          T0=296.15)
          annotation (Placement(transformation(extent={{-52,-62},{-24,-30}})));
        Components.Controls.sens_PT sens_PT3(redeclare package medium =
              ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
              extent={{6,7},{-6,-7}},
              rotation=-90,
              origin={-61,-18})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT
          annotation (Placement(transformation(extent={{66,76},{46,96}})));
        Modelica.Blocks.Sources.Step step(
          offset=808.15,
          startTime=6000,
          height=-325)
          annotation (Placement(transformation(extent={{20,104},{32,116}})));
        Components.main_components.Tank_liquid tank_liquid(
          rou=800,
          Cp=2.34e3,
          Area=0.05,
          L0=0.5,
          T0=499.15)
          annotation (Placement(transformation(extent={{60,46},{80,66}})));
        Modelica.Blocks.Sources.Step step1(
          offset=0,
          height=325,
          startTime=6300)
          annotation (Placement(transformation(extent={{20,126},{32,138}})));
        Modelica.Blocks.Math.Add add
          annotation (Placement(transformation(extent={{50,120},{62,132}})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT1
          annotation (Placement(transformation(extent={{-96,34},{-76,54}})));
        Modelica.Blocks.Sources.Step step2(
          startTime=6000,
          height=-2.4,
          offset=366.15)
          annotation (Placement(transformation(extent={{-154,70},{-142,82}})));
        Modelica.Blocks.Sources.Step step3(
          offset=0,
          height=2.4,
          startTime=6250)
          annotation (Placement(transformation(extent={{-140,92},{-128,104}})));
        Modelica.Blocks.Math.Add add1
          annotation (Placement(transformation(extent={{-110,78},{-96,92}})));
        Components.source_sink_ports.OutFlow_liquid_useM outFlow_liquid_useM
          annotation (Placement(transformation(extent={{-64,94},{-44,114}})));
        Modelica.Blocks.Sources.Step step6(
          startTime=6000,
          offset=0.2981,
          height=-0.1870) annotation (Placement(transformation(extent={{-120,
                  114},{-108,126}})));
        Modelica.Blocks.Sources.Step step7(
          offset=0,
          height=0.1870,
          startTime=6300) annotation (Placement(transformation(extent={{-140,
                  138},{-128,150}})));
        Modelica.Blocks.Math.Add add2
          annotation (Placement(transformation(extent={{-80,120},{-70,130}})));
        Components.Controls.SuperHeat superHeat(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{58,8},{78,28}})));
      equation
        connect(Gas_oil_Hx.inflow1, Evaporater.OutFlow_liquid) annotation (Line(
              points={{-24.36,73.8},{-24.36,72.9},{-36,72.9},{-36,34.16},{-14,
                34.16}}, color={28,108,200}));
        connect(Evaporater.InFlow_liquid, pump_liquid.outFlow) annotation (Line(
              points={{17.68,34.44},{24.84,34.44},{24.84,36},{30,36}}, color={
                28,108,200}));
        connect(const.y, pump_liquid.u) annotation (Line(points={{26.5,51},{
                26.5,50.5},{39.9,50.5},{39.9,39.9}}, color={0,0,127}));
        connect(Evaporater.OutFlow_working, sens_PT.inFlow) annotation (Line(
              points={{17.68,23.52},{27.84,23.52},{27.84,18.2},{32,18.2}},
              color={0,0,255}));
        connect(sens_PT1.outFlow, Evaporater.InFlow_working) annotation (Line(
              points={{-26,24.2},{-20,24.2},{-20,23.24},{-14,23.24}}, color={0,
                0,255}));
        connect(preHeater.OutFlow_liquid, outFlow_liquid1.A) annotation (Line(
              points={{-67.84,-2},{-70,-2},{-70,-16},{-74,-16}}, color={28,108,
                200}));
        connect(preHeater.OutFlow_working, sens_PT1.inFlow) annotation (Line(
              points={{-59.48,17.8},{-59.74,17.8},{-59.74,24.2},{-46,24.2}},
              color={0,0,255}));
        connect(turbine_simple.OutFlow, Condenser.InFlow_working) annotation (
            Line(points={{96.44,-24.22},{96.95,-24.22},{96.95,-41.57},{40,
                -41.57}},
              color={0,0,255}));
        connect(outFlow_liquid2.A, Condenser.OutFlow_liquid) annotation (Line(
              points={{66,-64},{46,-64},{46,-52.88},{40,-52.88}}, color={28,108,
                200}));
        connect(sens_PT2.inFlow, Condenser.OutFlow_working) annotation (Line(
              points={{2,-41.8},{6,-41.8},{6,-41.86},{7.33,-41.86}}, color={0,0,
                255}));
        connect(inFlow_liquid2.A, Condenser.InFlow_liquid) annotation (Line(
              points={{12.6,-80.2},{-7.3,-80.2},{-7.3,-53.17},{7.33,-53.17}},
              color={28,108,200}));
        connect(tank.InFlow, sens_PT2.outFlow) annotation (Line(points={{-29.6,
                -39.6},{-30,-39.6},{-30,-41.8},{-18,-41.8}}, color={0,0,255}));
        connect(tank.OutFlow, pump_type2_1.InFlow) annotation (Line(points={{
                -45.56,-54},{-45.56,-53.5},{-63.54,-53.5},{-63.54,-47.36}},
              color={0,0,255}));
        connect(sens_PT3.outFlow, preHeater.InFlow_working) annotation (Line(
              points={{-61.14,-12},{-59.26,-12},{-59.26,-2}}, color={0,0,255}));
        connect(pump_type2_1.OutFlow, sens_PT3.inFlow) annotation (Line(points=
                {{-73.16,-37.56},{-73.16,-32.7},{-61.14,-32.7},{-61.14,-24}},
              color={0,0,255}));
        connect(inFlow_liquid_useT.A, Gas_oil_Hx.inflow2) annotation (Line(
              points={{46.6,85.8},{26.3,85.8},{26.3,85.24},{11.64,85.24}},
              color={28,108,200}));
        connect(Gas_oil_Hx.outflow1, tank_liquid.inFlow) annotation (Line(
              points={{11.64,73.8},{80.82,73.8},{80.82,60.2},{76.2,60.2}},
              color={28,108,200}));
        connect(tank_liquid.outFlow, pump_liquid.inFlow) annotation (Line(
              points={{65.2,50.4},{65.2,35.2},{50,35.2},{50,36}}, color={28,108,
                200}));
        connect(step1.y, add.u1) annotation (Line(points={{32.6,132},{40,132},{
                40,129.6},{48.8,129.6}}, color={0,0,127}));
        connect(step.y, add.u2) annotation (Line(points={{32.6,110},{40,110},{
                40,122.4},{48.8,122.4}}, color={0,0,127}));
        connect(inFlow_liquid_useT1.A, preHeater.InFlow_liquid) annotation (
            Line(points={{-76.6,43.8},{-76.6,30.9},{-68.06,30.9},{-68.06,17.8}},
              color={28,108,200}));
        connect(step3.y, add1.u1) annotation (Line(points={{-127.4,98},{-120,98},
                {-120,89.2},{-111.4,89.2}}, color={0,0,127}));
        connect(step2.y, add1.u2) annotation (Line(points={{-141.4,76},{-120,76},
                {-120,80.8},{-111.4,80.8}}, color={0,0,127}));
        connect(add.y, inFlow_liquid_useT.u) annotation (Line(points={{62.6,126},
                {62.6,110.5},{57.2,110.5},{57.2,90}}, color={0,0,127}));
        connect(add1.y, inFlow_liquid_useT1.u) annotation (Line(points={{-95.3,
                85},{-95.3,66.5},{-87.2,66.5},{-87.2,48}}, color={0,0,127}));
        connect(step7.y, add2.u1) annotation (Line(points={{-127.4,144},{-96,
                144},{-96,128},{-81,128}}, color={0,0,127}));
        connect(step6.y, add2.u2) annotation (Line(points={{-107.4,120},{-96,
                120},{-96,122},{-81,122}}, color={0,0,127}));
        connect(outFlow_liquid_useM.A, Gas_oil_Hx.outflow2) annotation (Line(
              points={{-64,104},{-68,104},{-68,85.5},{-24,85.5}}, color={28,108,
                200}));
        connect(add2.y, outFlow_liquid_useM.u) annotation (Line(points={{-69.5,
                125},{-54.65,125},{-54.65,108.2},{-54.2,108.2}}, color={0,0,127}));
        connect(sens_PT.outFlow, superHeat.inFlow) annotation (Line(points={{52,
                18.2},{56,18.2},{58,18.2}}, color={0,0,255}));
        connect(superHeat.outFlow, turbine_simple.InFlow) annotation (Line(
              points={{78,18.2},{82,18.2},{82,-5.86},{83,-5.86}}, color={0,0,
                255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-48,10},{48,-34}},
                lineColor={28,108,200},
                textString="内燃机余热回收")}));
      end CycleOneICE_part1;

      model CycleOneICE_part2
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end CycleOneICE_part2;

      model CycleICE_five_condition
        Components.main_components.Hx1counter_liquid Gas_oil_Hx(
          Cp1=2.34e3,
          U1=1300,
          U2=800,
          Cp_wall=400,
          M_wall=20,
          F_liquid1=3,
          F_liquid2=3,
          V_liquid1=0.7,
          V_liquid2=0.7,
          M1=500,
          N=20,
          T1_liquid(displayUnit="degC") = 447.15,
          T1n_liquid(displayUnit="degC") = 499.55,
          T2_liquid(displayUnit="degC") = 808.15,
          T2n_liquid(displayUnit="degC") = 455.42,
          Cp2=1.109e3,
          M2=0.5)
          annotation (Placement(transformation(extent={{-24,54},{12,80}})));
        Components.main_components.Hx1counter Evaporater(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Utp=3000,
          Uv=800,
          U=1000,
          M=500,
          Cp=2.34e3,
          M_wall=20,
          Cp_wall=400,
          F_working=2.3,
          F_liquid=2.3,
          m_liquid=1,
          N=20,
          m_working=0.7,
          Tin_start=365.33,
          Tout_start=393.15,
          T1_liquid=499.55,
          Tn_liquid=447.85,
          pstart=1895000)
          annotation (Placement(transformation(extent={{-14,2},{18,30}})));
        Components.main_components.Pump_liquid pump_liquid
          annotation (Placement(transformation(extent={{50,14},{30,34}})));
        Modelica.Blocks.Sources.Constant const(k=1)
          annotation (Placement(transformation(extent={{16,34},{26,44}})));
        Components.Controls.sens_PT sens_B(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{32,-4},{52,16}})));
        Components.Controls.sens_PT sens_C(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{-46,2},{-26,22}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(m0=2.1)
          annotation (Placement(transformation(extent={{-78,-38},{-98,-18}})));
        Components.main_components.Hx1counter preHeater(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Utp=3000,
          Uv=800,
          U=1000,
          M_wall=20,
          Cp_wall=400,
          N=20,
          M=700,
          Cp=4.2e3,
          m_liquid=2.1,
          F_working=10,
          F_liquid=10,
          m_working=0.7,
          Tin_start=296.45,
          Tout_start=365.35,
          T1_liquid=366.15,
          Tn_liquid=358.45,
          pstart=1895000) annotation (Placement(transformation(
              extent={{-10,-11},{10,11}},
              rotation=90,
              origin={-63,-4})));
        Components.main_components.Turbine_simple turbine_simple(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          Np=35,
          V_s=1.904e-4,
          pin_start=1895000,
          pout_start=176000,
          Tin_start=398.15)
          annotation (Placement(transformation(extent={{76,-42},{104,-8}})));
        Components.main_components.Hx1counter Condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Utp=3000,
          Uv=800,
          U=1000,
          M_wall=20,
          Cp_wall=400,
          N=20,
          M=700,
          Cp=4.2e3,
          m_working=0.7,
          m_liquid=10,
          F_working=20,
          F_liquid=20,
          Tin_start=347.11,
          Tout_start=295.38,
          T1_liquid=290.15,
          Tn_liquid=294.15,
          pstart=176000) annotation (Placement(transformation(
              extent={{-16.5,-14.5},{16.5,14.5}},
              rotation=180,
              origin={23.5,-58.5})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid2(m0=15)
          annotation (Placement(transformation(extent={{66,-86},{86,-66}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid2(T0(
              displayUnit="degC") = 290.15)
          annotation (Placement(transformation(extent={{32,-102},{12,-82}})));
        Components.Controls.sens_PT sens_PT2(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{2,-64},{-18,-44}})));
        Components.main_components.Tank tank(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          Aera=0.065,
          L0=0.65,
          p0=176000,
          T0=296.15)
          annotation (Placement(transformation(extent={{-52,-74},{-24,-42}})));
        Components.Controls.sens_PT sens_D(redeclare package medium =
              ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
              extent={{6,7},{-6,-7}},
              rotation=-90,
              origin={-61,-30})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT
          annotation (Placement(transformation(extent={{70,64},{50,84}})));
        Components.main_components.Tank_liquid tank_liquid(
          rou=800,
          Cp=2.34e3,
          Area=0.05,
          L0=0.5,
          T0=499.15)
          annotation (Placement(transformation(extent={{60,34},{80,54}})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT1
          annotation (Placement(transformation(extent={{-96,22},{-76,42}})));
        Components.source_sink_ports.OutFlow_liquid_useM outFlow_liquid_useM
          annotation (Placement(transformation(extent={{-64,82},{-44,102}})));
        Components.Controls.SuperHeat superHeat(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{58,-4},{78,16}})));
        Components.main_components.Pump_type4 pump_type4_1(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          v_s=1.0395e-5,
          p0=176000,
          pstart_out=1895000,
          T0=295.15,
          v_max=1.9e-5)
          annotation (Placement(transformation(extent={{-52,-76},{-88,-40}})));
        Modelica.Blocks.Sources.Constant const2(k=50)
          annotation (Placement(transformation(extent={{-102,-62},{-92,-52}})));
        Modelica.Blocks.Sources.Constant const1(k=365.65)
          annotation (Placement(transformation(extent={{-106,46},{-96,56}})));
        Modelica.Blocks.Sources.Step step1_2(
          height=-0.0229,
          offset=0.2981,
          startTime=8000) annotation (Placement(transformation(extent={{-122,
                  104},{-102,124}})));
        Modelica.Blocks.Sources.Step step1_3(
          offset=0.2981,
          height=-0.0395,
          startTime=8000) annotation (Placement(transformation(extent={{-124,
                  136},{-104,156}})));
        Modelica.Blocks.Sources.Step step1_4(
          offset=0.2981,
          height=-0.0746,
          startTime=8000)
          annotation (Placement(transformation(extent={{-132,172},{-112,192}})));
        Modelica.Blocks.Sources.Step step1_5(
          offset=0.2981,
          height=-0.1284,
          startTime=8000) annotation (Placement(transformation(extent={{-120,
                  208},{-100,228}})));
        Modelica.Blocks.Sources.Step stepT1_2(
          offset=808.15,
          height=-16,
          startTime=8000)
          annotation (Placement(transformation(extent={{20,96},{40,116}})));
        Modelica.Blocks.Sources.Step stepT1_3(
          offset=808.15,
          height=-37,
          startTime=8000)
          annotation (Placement(transformation(extent={{16,130},{36,150}})));
        Modelica.Blocks.Sources.Step stepT1_4(
          offset=808.15,
          height=-61,
          startTime=8000)
          annotation (Placement(transformation(extent={{20,166},{40,186}})));
        Modelica.Blocks.Sources.Step stepT1_5(
          offset=808.15,
          height=-115,
          startTime=8000)
          annotation (Placement(transformation(extent={{18,202},{38,222}})));
      equation
        connect(Gas_oil_Hx.inflow1,Evaporater. OutFlow_liquid) annotation (Line(
              points={{-24.36,61.8},{-24.36,60.9},{-36,60.9},{-36,22.16},{-14,
                22.16}}, color={28,108,200}));
        connect(Evaporater.InFlow_liquid,pump_liquid. outFlow) annotation (Line(
              points={{17.68,22.44},{24.84,22.44},{24.84,24},{30,24}}, color={
                28,108,200}));
        connect(const.y,pump_liquid. u) annotation (Line(points={{26.5,39},{
                26.5,38.5},{39.9,38.5},{39.9,27.9}}, color={0,0,127}));
        connect(Evaporater.OutFlow_working, sens_B.inFlow) annotation (Line(
              points={{17.68,11.52},{27.84,11.52},{27.84,6.2},{32,6.2}}, color=
                {0,0,255}));
        connect(sens_C.outFlow, Evaporater.InFlow_working) annotation (Line(
              points={{-26,12.2},{-20,12.2},{-20,11.24},{-14,11.24}}, color={0,
                0,255}));
        connect(preHeater.OutFlow_liquid,outFlow_liquid1. A) annotation (Line(
              points={{-67.84,-14},{-70,-14},{-70,-28},{-78,-28}},
                                                                 color={28,108,
                200}));
        connect(preHeater.OutFlow_working, sens_C.inFlow) annotation (Line(
              points={{-59.48,5.8},{-59.74,5.8},{-59.74,12.2},{-46,12.2}},
              color={0,0,255}));
        connect(turbine_simple.OutFlow,Condenser. InFlow_working) annotation (
            Line(points={{96.44,-36.22},{96.95,-36.22},{96.95,-53.57},{40,
                -53.57}},
              color={0,0,255}));
        connect(outFlow_liquid2.A,Condenser. OutFlow_liquid) annotation (Line(
              points={{66,-76},{46,-76},{46,-64.88},{40,-64.88}}, color={28,108,
                200}));
        connect(sens_PT2.inFlow,Condenser. OutFlow_working) annotation (Line(
              points={{2,-53.8},{6,-53.8},{6,-53.86},{7.33,-53.86}}, color={0,0,
                255}));
        connect(inFlow_liquid2.A,Condenser. InFlow_liquid) annotation (Line(
              points={{12.6,-92.2},{-7.3,-92.2},{-7.3,-65.17},{7.33,-65.17}},
              color={28,108,200}));
        connect(tank.InFlow,sens_PT2. outFlow) annotation (Line(points={{-29.6,
                -51.6},{-30,-51.6},{-30,-53.8},{-18,-53.8}}, color={0,0,255}));
        connect(sens_D.outFlow, preHeater.InFlow_working) annotation (Line(
              points={{-61.14,-24},{-59.26,-24},{-59.26,-14}}, color={0,0,255}));
        connect(inFlow_liquid_useT.A,Gas_oil_Hx. inflow2) annotation (Line(
              points={{50.6,73.8},{26.3,73.8},{26.3,73.24},{11.64,73.24}},
              color={28,108,200}));
        connect(Gas_oil_Hx.outflow1,tank_liquid. inFlow) annotation (Line(
              points={{11.64,61.8},{80.82,61.8},{80.82,48.2},{76.2,48.2}},
              color={28,108,200}));
        connect(tank_liquid.outFlow,pump_liquid. inFlow) annotation (Line(
              points={{65.2,38.4},{65.2,23.2},{50,23.2},{50,24}}, color={28,108,
                200}));
        connect(inFlow_liquid_useT1.A,preHeater. InFlow_liquid) annotation (
            Line(points={{-76.6,31.8},{-76.6,18.9},{-68.06,18.9},{-68.06,5.8}},
              color={28,108,200}));
        connect(outFlow_liquid_useM.A,Gas_oil_Hx. outflow2) annotation (Line(
              points={{-64,92},{-68,92},{-68,73.5},{-24,73.5}},   color={28,108,
                200}));
        connect(sens_B.outFlow, superHeat.inFlow) annotation (Line(points={{52,
                6.2},{56,6.2},{58,6.2}}, color={0,0,255}));
        connect(superHeat.outFlow, turbine_simple.InFlow) annotation (Line(
              points={{78,6.2},{82,6.2},{82,-17.86},{83,-17.86}}, color={0,0,
                255}));
        connect(pump_type4_1.InFlow, tank.OutFlow) annotation (Line(points={{
                -62.44,-62.32},{-53.22,-62.32},{-53.22,-66},{-45.56,-66}},
              color={0,0,255}));
        connect(pump_type4_1.OutFlow, sens_D.inFlow) annotation (Line(points={{
                -75.76,-49.72},{-67.88,-49.72},{-67.88,-36},{-61.14,-36}},
              color={0,0,255}));
        connect(const2.y, pump_type4_1.u) annotation (Line(points={{-91.5,-57},
                {-82.75,-57},{-82.75,-55.84},{-76.12,-55.84}}, color={0,0,127}));
        connect(const1.y, inFlow_liquid_useT1.u) annotation (Line(points={{
                -95.5,51},{-95.5,51.5},{-87.2,51.5},{-87.2,36}}, color={0,0,127}));
        connect(step1_5.y, outFlow_liquid_useM.u) annotation (Line(points={{-99,
                218},{-76,218},{-76,96.2},{-54.2,96.2}}, color={0,0,127}));
        connect(stepT1_5.y, inFlow_liquid_useT.u) annotation (Line(points={{39,
                212},{50,212},{50,78},{61.2,78}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-48,0},{48,-44}},
                lineColor={28,108,200},
                textString="内燃机余热回收")}));
      end CycleICE_five_condition;

      model Cycle_five_condition_ORC
        Components.main_components.Hx1counter Evaporater(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Utp=3000,
          Uv=800,
          U=1000,
          M_wall=20,
          Cp_wall=400,
          F_working=2.3,
          F_liquid=2.3,
          m_liquid=1,
          N=20,
          m_working=0.7,
          Cp=1.109e3,
          M=0.5,
          Tin_start=365.33,
          Tout_start=423.15,
          T1_liquid=808.15,
          Tn_liquid=382.15,
          pstart=1895000)
          annotation (Placement(transformation(extent={{-4,22},{28,50}})));
        Components.Controls.sens_PT sens_B(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{40,-4},{60,16}})));
        Components.Controls.sens_PT sens_C(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{-38,2},{-18,22}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(m0=2.1)
          annotation (Placement(transformation(extent={{-70,-38},{-90,-18}})));
        Components.main_components.Hx1counter preHeater(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Utp=3000,
          Uv=800,
          U=1000,
          M_wall=20,
          Cp_wall=400,
          N=20,
          M=700,
          Cp=4.2e3,
          m_liquid=2.1,
          F_working=10,
          F_liquid=10,
          m_working=0.7,
          Tin_start=296.45,
          Tout_start=365.35,
          T1_liquid=366.15,
          Tn_liquid=358.45,
          pstart=1895000) annotation (Placement(transformation(
              extent={{-10,-11},{10,11}},
              rotation=90,
              origin={-55,-4})));
        Components.main_components.Turbine_simple turbine_simple(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          V_s=1.904e-4,
          Np=35,
          pin_start=1895000,
          pout_start=176000,
          Tin_start=398.15)
          annotation (Placement(transformation(extent={{84,-44},{112,-10}})));
        Components.main_components.Hx1counter Condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Utp=3000,
          Uv=800,
          U=1000,
          M_wall=20,
          Cp_wall=400,
          N=20,
          M=700,
          Cp=4.2e3,
          m_working=0.7,
          m_liquid=10,
          F_working=20,
          F_liquid=20,
          Tin_start=347.11,
          Tout_start=295.38,
          T1_liquid=290.15,
          Tn_liquid=294.15,
          pstart=176000) annotation (Placement(transformation(
              extent={{-16.5,-14.5},{16.5,14.5}},
              rotation=180,
              origin={31.5,-58.5})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid2(m0=22)
          annotation (Placement(transformation(extent={{74,-86},{94,-66}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid2(T0(
              displayUnit="degC") = 289.15)
          annotation (Placement(transformation(extent={{40,-102},{20,-82}})));
        Components.Controls.sens_PT sens_PT2(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{10,-64},{-10,-44}})));
        Components.Controls.sens_PT sens_D(redeclare package medium =
              ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
              extent={{6,7},{-6,-7}},
              rotation=-90,
              origin={-53,-30})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT1
          annotation (Placement(transformation(extent={{-88,22},{-68,42}})));
        Components.Controls.SuperHeat superHeat(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{66,-4},{86,16}})));
        Components.main_components.Pump_type4 pump_type4_1(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          v_s=1.0395e-5,
          p0=176000,
          pstart_out=1895000,
          T0=295.15,
          v_max=1.9e-5)
          annotation (Placement(transformation(extent={{-44,-76},{-80,-40}})));
        Modelica.Blocks.Sources.Constant const2(k=50)
          annotation (Placement(transformation(extent={{-94,-62},{-84,-52}})));
        Modelica.Blocks.Sources.Constant const1(k=365.65)
          annotation (Placement(transformation(extent={{-98,46},{-88,56}})));
        Components.source_sink_ports.OutFlow_liquid_useM outFlow_liquid_useM
          annotation (Placement(transformation(extent={{-28,44},{-48,64}})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT
          annotation (Placement(transformation(extent={{72,32},{52,52}})));
        Modelica.Blocks.Sources.Step step1_2(
          height=-0.0229,
          offset=0.2981,
          startTime=4000)
          annotation (Placement(transformation(extent={{-76,82},{-56,102}})));
        Modelica.Blocks.Sources.Step stepT1_2(
          offset=808.15,
          height=-16,
          startTime=4000)
          annotation (Placement(transformation(extent={{32,84},{52,104}})));
        Components.Controls.SuperHeat superHeat1(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{-22,-90},{-42,-70}})));
        Components.main_components.Tank_rou tank_rou(
          L0=0.65,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0=176000,
          T0=296.15,
          Aera=0.5)
          annotation (Placement(transformation(extent={{-66,-94},{-46,-74}})));
        Modelica.Blocks.Sources.Step step1_3(
          offset=0.2981,
          height=-0.0395,
          startTime=4000)
          annotation (Placement(transformation(extent={{-88,120},{-68,140}})));
        Modelica.Blocks.Sources.Step step1_4(
          offset=0.2981,
          height=-0.0746,
          startTime=4000)
          annotation (Placement(transformation(extent={{-76,156},{-56,176}})));
        Modelica.Blocks.Sources.Step step1_5(
          offset=0.2981,
          height=-0.1284,
          startTime=4000)
          annotation (Placement(transformation(extent={{-72,188},{-52,208}})));
        Modelica.Blocks.Sources.Step stepT1_3(
          offset=808.15,
          height=-37,
          startTime=4000)
          annotation (Placement(transformation(extent={{34,116},{54,136}})));
        Modelica.Blocks.Sources.Step stepT1_4(
          offset=808.15,
          height=-61,
          startTime=4000)
          annotation (Placement(transformation(extent={{32,150},{52,170}})));
        Modelica.Blocks.Sources.Step stepT1_5(
          offset=808.15,
          height=-115,
          startTime=4000)
          annotation (Placement(transformation(extent={{32,190},{52,210}})));
      equation
        connect(Evaporater.OutFlow_working, sens_B.inFlow) annotation (Line(
              points={{27.68,31.52},{35.84,31.52},{35.84,6.2},{40,6.2}}, color=
                {0,0,255}));
        connect(sens_C.outFlow, Evaporater.InFlow_working) annotation (Line(
              points={{-18,12.2},{-12,12.2},{-12,31.24},{-4,31.24}}, color={0,0,
                255}));
        connect(preHeater.OutFlow_liquid,outFlow_liquid1. A) annotation (Line(
              points={{-59.84,-14},{-62,-14},{-62,-28},{-70,-28}},
                                                                 color={28,108,
                200}));
        connect(preHeater.OutFlow_working, sens_C.inFlow) annotation (Line(
              points={{-51.48,5.8},{-51.74,5.8},{-51.74,12.2},{-38,12.2}},
              color={0,0,255}));
        connect(turbine_simple.OutFlow,Condenser. InFlow_working) annotation (
            Line(points={{104.44,-38.22},{104.95,-38.22},{104.95,-53.57},{48,
                -53.57}},
              color={0,0,255}));
        connect(outFlow_liquid2.A,Condenser. OutFlow_liquid) annotation (Line(
              points={{74,-76},{54,-76},{54,-64.88},{48,-64.88}}, color={28,108,
                200}));
        connect(sens_PT2.inFlow,Condenser. OutFlow_working) annotation (Line(
              points={{10,-53.8},{14,-53.8},{14,-53.86},{15.33,-53.86}},
                                                                     color={0,0,
                255}));
        connect(inFlow_liquid2.A,Condenser. InFlow_liquid) annotation (Line(
              points={{20.6,-92.2},{0.7,-92.2},{0.7,-65.17},{15.33,-65.17}},
              color={28,108,200}));
        connect(sens_D.outFlow, preHeater.InFlow_working) annotation (Line(
              points={{-53.14,-24},{-51.26,-24},{-51.26,-14}}, color={0,0,255}));
        connect(inFlow_liquid_useT1.A,preHeater. InFlow_liquid) annotation (
            Line(points={{-68.6,31.8},{-68.6,18.9},{-60.06,18.9},{-60.06,5.8}},
              color={28,108,200}));
        connect(sens_B.outFlow, superHeat.inFlow) annotation (Line(points={{60,
                6.2},{64,6.2},{66,6.2}}, color={0,0,255}));
        connect(superHeat.outFlow, turbine_simple.InFlow) annotation (Line(
              points={{86,6.2},{90,6.2},{90,-19.86},{91,-19.86}}, color={0,0,
                255}));
        connect(pump_type4_1.OutFlow, sens_D.inFlow) annotation (Line(points={{
                -67.76,-49.72},{-59.88,-49.72},{-59.88,-36},{-53.14,-36}},
              color={0,0,255}));
        connect(const2.y, pump_type4_1.u) annotation (Line(points={{-83.5,-57},
                {-74.75,-57},{-74.75,-55.84},{-68.12,-55.84}}, color={0,0,127}));
        connect(const1.y, inFlow_liquid_useT1.u) annotation (Line(points={{
                -87.5,51},{-87.5,51.5},{-79.2,51.5},{-79.2,36}}, color={0,0,127}));
        connect(inFlow_liquid_useT.A, Evaporater.InFlow_liquid) annotation (
            Line(points={{52.6,41.8},{38.3,41.8},{38.3,42.44},{27.68,42.44}},
              color={28,108,200}));
        connect(outFlow_liquid_useM.A, Evaporater.OutFlow_liquid) annotation (
            Line(points={{-28,54},{-18,54},{-18,42.16},{-4,42.16}}, color={28,
                108,200}));
        connect(sens_PT2.outFlow, superHeat1.inFlow) annotation (Line(points={{
                -10,-53.8},{-16,-53.8},{-16,-79.8},{-22,-79.8}}, color={0,0,255}));
        connect(tank_rou.InFlow, superHeat1.outFlow) annotation (Line(points={{
                -50,-80},{-42,-80},{-42,-79.8}}, color={0,0,255}));
        connect(pump_type4_1.InFlow, tank_rou.OutFlow) annotation (Line(points=
                {{-54.44,-62.32},{-54.44,-76.16},{-61.4,-76.16},{-61.4,-89}},
              color={0,0,255}));
        connect(step1_5.y, outFlow_liquid_useM.u) annotation (Line(points={{-51,
                198},{-44,198},{-44,58.2},{-37.8,58.2}}, color={0,0,127}));
        connect(stepT1_5.y, inFlow_liquid_useT.u) annotation (Line(points={{53,
                200},{58,200},{58,46},{63.2,46}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-40,0},{56,-44}},
                lineColor={28,108,200},
                textString="内燃机余热回收")}));
      end Cycle_five_condition_ORC;

      model Cycle_condition_change_OS_ORC
        Components.main_components.Hx1counter_liquid Gas_oil_Hx(
          Cp1=2.34e3,
          U1=1300,
          U2=800,
          Cp_wall=400,
          M_wall=20,
          F_liquid1=3,
          F_liquid2=3,
          V_liquid1=0.7,
          V_liquid2=0.7,
          M1=500,
          N=20,
          T1_liquid(displayUnit="degC") = 447.15,
          T1n_liquid(displayUnit="degC") = 499.55,
          T2_liquid(displayUnit="degC") = 808.15,
          T2n_liquid(displayUnit="degC") = 455.42,
          Cp2=1.109e3,
          M2=0.5)
          annotation (Placement(transformation(extent={{-24,48},{12,74}})));
        Components.main_components.Hx1counter Evaporater(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Utp=3000,
          Uv=800,
          U=1000,
          M=500,
          Cp=2.34e3,
          M_wall=20,
          Cp_wall=400,
          F_working=2.3,
          F_liquid=2.3,
          m_liquid=1,
          N=20,
          m_working=0.7,
          Tin_start=365.33,
          Tout_start=393.15,
          T1_liquid=499.55,
          Tn_liquid=447.85,
          pstart=1895000)
          annotation (Placement(transformation(extent={{-14,-4},{18,24}})));
        Components.main_components.Pump_liquid pump_liquid
          annotation (Placement(transformation(extent={{50,8},{30,28}})));
        Modelica.Blocks.Sources.Constant const(k=1)
          annotation (Placement(transformation(extent={{16,28},{26,38}})));
        Components.Controls.sens_PT sens_B(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{32,-10},{52,10}})));
        Components.Controls.sens_PT sens_C(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{-46,-4},{-26,16}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(m0=2.1)
          annotation (Placement(transformation(extent={{-78,-44},{-98,-24}})));
        Components.main_components.Hx1counter preHeater(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Utp=3000,
          Uv=800,
          U=1000,
          M_wall=20,
          Cp_wall=400,
          N=20,
          M=700,
          Cp=4.2e3,
          m_liquid=2.1,
          F_working=10,
          F_liquid=10,
          m_working=0.7,
          Tin_start=296.45,
          Tout_start=365.35,
          T1_liquid=366.15,
          Tn_liquid=358.45,
          pstart=1895000) annotation (Placement(transformation(
              extent={{-10,-11},{10,11}},
              rotation=90,
              origin={-63,-10})));
        Components.main_components.Turbine_simple turbine_simple(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          Np=35,
          V_s=1.904e-4,
          pin_start=1895000,
          pout_start=176000,
          Tin_start=398.15)
          annotation (Placement(transformation(extent={{76,-48},{104,-14}})));
        Components.main_components.Hx1counter Condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Utp=3000,
          Uv=800,
          U=1000,
          M_wall=20,
          Cp_wall=400,
          N=20,
          M=700,
          Cp=4.2e3,
          m_working=0.7,
          m_liquid=10,
          F_working=20,
          F_liquid=20,
          Tin_start=347.11,
          Tout_start=295.38,
          T1_liquid=290.15,
          Tn_liquid=294.15,
          pstart=176000) annotation (Placement(transformation(
              extent={{-16.5,-14.5},{16.5,14.5}},
              rotation=180,
              origin={23.5,-64.5})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid2(m0=15)
          annotation (Placement(transformation(extent={{66,-92},{86,-72}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid2(T0(
              displayUnit="degC") = 290.15)
          annotation (Placement(transformation(extent={{32,-108},{12,-88}})));
        Components.Controls.sens_PT sens_PT2(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{2,-70},{-18,-50}})));
        Components.main_components.Tank tank(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          Aera=0.065,
          L0=0.65,
          p0=176000,
          T0=296.15)
          annotation (Placement(transformation(extent={{-52,-80},{-24,-48}})));
        Components.Controls.sens_PT sens_D(redeclare package medium =
              ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
              extent={{6,7},{-6,-7}},
              rotation=-90,
              origin={-61,-36})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT
          annotation (Placement(transformation(extent={{70,58},{50,78}})));
        Components.main_components.Tank_liquid tank_liquid(
          rou=800,
          Cp=2.34e3,
          Area=0.05,
          L0=0.5,
          T0=499.15)
          annotation (Placement(transformation(extent={{60,28},{80,48}})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT1
          annotation (Placement(transformation(extent={{-96,16},{-76,36}})));
        Components.source_sink_ports.OutFlow_liquid_useM outFlow_liquid_useM
          annotation (Placement(transformation(extent={{-64,76},{-44,96}})));
        Components.Controls.SuperHeat superHeat(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{58,-10},{78,10}})));
        Components.main_components.Pump_type4 pump_type4_1(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          v_s=1.0395e-5,
          p0=176000,
          pstart_out=1895000,
          T0=295.15,
          v_max=1.9e-5)
          annotation (Placement(transformation(extent={{-52,-82},{-88,-46}})));
        Modelica.Blocks.Sources.Constant const2(k=50)
          annotation (Placement(transformation(extent={{-102,-68},{-92,-58}})));
        Modelica.Blocks.Sources.Constant const1(k=365.65)
          annotation (Placement(transformation(extent={{-106,40},{-96,50}})));
        Modelica.Blocks.Sources.Step step1_5(
          offset=0.2981,
          height=-0.1284,
          startTime=8000) annotation (Placement(transformation(extent={{-132,
                  140},{-112,160}})));
        Modelica.Blocks.Sources.Step stepT1_5(
          offset=808.15,
          height=-115,
          startTime=8000)
          annotation (Placement(transformation(extent={{8,144},{28,164}})));
        Modelica.Blocks.Sources.Step step(height=0.1284, startTime=8600)
          annotation (Placement(transformation(extent={{-132,106},{-112,126}})));
        Modelica.Blocks.Sources.Step stepT1_1(
          height=115,
          offset=0,
          startTime=8600)
          annotation (Placement(transformation(extent={{6,100},{26,120}})));
        Modelica.Blocks.Math.Add add
          annotation (Placement(transformation(extent={{-90,118},{-70,138}})));
        Modelica.Blocks.Math.Add add1
          annotation (Placement(transformation(extent={{40,124},{60,144}})));
        Modelica.Blocks.Sources.Constant const3(k=0.2981)
          annotation (Placement(transformation(extent={{-92,72},{-82,82}})));
        Modelica.Blocks.Sources.Constant const4(k=808.15)
          annotation (Placement(transformation(extent={{86,68},{96,78}})));
      equation
        connect(Gas_oil_Hx.inflow1,Evaporater. OutFlow_liquid) annotation (Line(
              points={{-24.36,55.8},{-24.36,54.9},{-36,54.9},{-36,16.16},{-14,
                16.16}}, color={28,108,200}));
        connect(Evaporater.InFlow_liquid,pump_liquid. outFlow) annotation (Line(
              points={{17.68,16.44},{24.84,16.44},{24.84,18},{30,18}}, color={
                28,108,200}));
        connect(const.y,pump_liquid. u) annotation (Line(points={{26.5,33},{
                26.5,32.5},{39.9,32.5},{39.9,21.9}}, color={0,0,127}));
        connect(Evaporater.OutFlow_working, sens_B.inFlow) annotation (Line(
              points={{17.68,5.52},{27.84,5.52},{27.84,0.2},{32,0.2}}, color={0,
                0,255}));
        connect(sens_C.outFlow, Evaporater.InFlow_working) annotation (Line(
              points={{-26,6.2},{-20,6.2},{-20,5.24},{-14,5.24}}, color={0,0,
                255}));
        connect(preHeater.OutFlow_liquid,outFlow_liquid1. A) annotation (Line(
              points={{-67.84,-20},{-70,-20},{-70,-34},{-78,-34}},
                                                                 color={28,108,
                200}));
        connect(preHeater.OutFlow_working, sens_C.inFlow) annotation (Line(
              points={{-59.48,-0.2},{-59.74,-0.2},{-59.74,6.2},{-46,6.2}},
              color={0,0,255}));
        connect(turbine_simple.OutFlow,Condenser. InFlow_working) annotation (
            Line(points={{96.44,-42.22},{96.95,-42.22},{96.95,-59.57},{40,
                -59.57}},
              color={0,0,255}));
        connect(outFlow_liquid2.A,Condenser. OutFlow_liquid) annotation (Line(
              points={{66,-82},{46,-82},{46,-70.88},{40,-70.88}}, color={28,108,
                200}));
        connect(sens_PT2.inFlow,Condenser. OutFlow_working) annotation (Line(
              points={{2,-59.8},{6,-59.8},{6,-59.86},{7.33,-59.86}}, color={0,0,
                255}));
        connect(inFlow_liquid2.A,Condenser. InFlow_liquid) annotation (Line(
              points={{12.6,-98.2},{-7.3,-98.2},{-7.3,-71.17},{7.33,-71.17}},
              color={28,108,200}));
        connect(tank.InFlow,sens_PT2. outFlow) annotation (Line(points={{-29.6,
                -57.6},{-30,-57.6},{-30,-59.8},{-18,-59.8}}, color={0,0,255}));
        connect(sens_D.outFlow, preHeater.InFlow_working) annotation (Line(
              points={{-61.14,-30},{-59.26,-30},{-59.26,-20}}, color={0,0,255}));
        connect(inFlow_liquid_useT.A,Gas_oil_Hx. inflow2) annotation (Line(
              points={{50.6,67.8},{26.3,67.8},{26.3,67.24},{11.64,67.24}},
              color={28,108,200}));
        connect(Gas_oil_Hx.outflow1,tank_liquid. inFlow) annotation (Line(
              points={{11.64,55.8},{80.82,55.8},{80.82,42.2},{76.2,42.2}},
              color={28,108,200}));
        connect(tank_liquid.outFlow,pump_liquid. inFlow) annotation (Line(
              points={{65.2,32.4},{65.2,17.2},{50,17.2},{50,18}}, color={28,108,
                200}));
        connect(inFlow_liquid_useT1.A,preHeater. InFlow_liquid) annotation (
            Line(points={{-76.6,25.8},{-76.6,12.9},{-68.06,12.9},{-68.06,-0.2}},
              color={28,108,200}));
        connect(outFlow_liquid_useM.A,Gas_oil_Hx. outflow2) annotation (Line(
              points={{-64,86},{-68,86},{-68,67.5},{-24,67.5}},   color={28,108,
                200}));
        connect(sens_B.outFlow, superHeat.inFlow) annotation (Line(points={{52,
                0.2},{56,0.2},{58,0.2}}, color={0,0,255}));
        connect(superHeat.outFlow, turbine_simple.InFlow) annotation (Line(
              points={{78,0.2},{82,0.2},{82,-23.86},{83,-23.86}}, color={0,0,
                255}));
        connect(pump_type4_1.InFlow, tank.OutFlow) annotation (Line(points={{
                -62.44,-68.32},{-53.22,-68.32},{-53.22,-72},{-45.56,-72}},
              color={0,0,255}));
        connect(pump_type4_1.OutFlow, sens_D.inFlow) annotation (Line(points={{
                -75.76,-55.72},{-67.88,-55.72},{-67.88,-42},{-61.14,-42}},
              color={0,0,255}));
        connect(const2.y, pump_type4_1.u) annotation (Line(points={{-91.5,-63},
                {-82.75,-63},{-82.75,-61.84},{-76.12,-61.84}}, color={0,0,127}));
        connect(const1.y, inFlow_liquid_useT1.u) annotation (Line(points={{
                -95.5,45},{-95.5,45.5},{-87.2,45.5},{-87.2,30}}, color={0,0,127}));
        connect(step1_5.y, add.u1) annotation (Line(points={{-111,150},{-100.5,
                150},{-100.5,134},{-92,134}}, color={0,0,127}));
        connect(step.y, add.u2) annotation (Line(points={{-111,116},{-102,116},
                {-102,122},{-92,122}}, color={0,0,127}));
        connect(stepT1_5.y, add1.u1) annotation (Line(points={{29,154},{36,154},
                {36,140},{38,140}}, color={0,0,127}));
        connect(stepT1_1.y, add1.u2) annotation (Line(points={{27,110},{36,110},
                {36,128},{38,128}}, color={0,0,127}));
        connect(add.y, outFlow_liquid_useM.u) annotation (Line(points={{-69,128},
                {-62,128},{-62,90.2},{-54.2,90.2}}, color={0,0,127}));
        connect(add1.y, inFlow_liquid_useT.u) annotation (Line(points={{61,134},
                {61.2,134},{61.2,72}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-48,-6},{48,-50}},
                lineColor={28,108,200},
                textString="内燃机余热回收")}));
      end Cycle_condition_change_OS_ORC;

      model Cycle_condition_change_ORC
        Components.main_components.Hx1counter Evaporater(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Utp=3000,
          Uv=800,
          U=1000,
          M_wall=20,
          Cp_wall=400,
          F_working=2.3,
          F_liquid=2.3,
          m_liquid=1,
          N=20,
          m_working=0.7,
          Cp=1.109e3,
          M=0.5,
          Tin_start=365.33,
          Tout_start=423.15,
          T1_liquid=808.15,
          Tn_liquid=382.15,
          pstart=1895000)
          annotation (Placement(transformation(extent={{-18,-8},{14,20}})));
        Components.Controls.sens_PT sens_B(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{26,-34},{46,-14}})));
        Components.Controls.sens_PT sens_C(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{-52,-28},{-32,-8}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(m0=2.1)
          annotation (Placement(transformation(extent={{-84,-68},{-104,-48}})));
        Components.main_components.Hx1counter preHeater(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Utp=3000,
          Uv=800,
          U=1000,
          M_wall=20,
          Cp_wall=400,
          N=20,
          M=700,
          Cp=4.2e3,
          m_liquid=2.1,
          F_working=10,
          F_liquid=10,
          m_working=0.7,
          Tin_start=296.45,
          Tout_start=365.35,
          T1_liquid=366.15,
          Tn_liquid=358.45,
          pstart=1895000) annotation (Placement(transformation(
              extent={{-10,-11},{10,11}},
              rotation=90,
              origin={-69,-34})));
        Components.main_components.Turbine_simple turbine_simple(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          V_s=1.904e-4,
          Np=35,
          pin_start=1895000,
          pout_start=176000,
          Tin_start=398.15)
          annotation (Placement(transformation(extent={{70,-72},{98,-38}})));
        Components.main_components.Hx1counter Condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Utp=3000,
          Uv=800,
          U=1000,
          M_wall=20,
          Cp_wall=400,
          N=20,
          M=700,
          Cp=4.2e3,
          m_working=0.7,
          m_liquid=10,
          F_working=20,
          F_liquid=20,
          Tin_start=347.11,
          Tout_start=295.38,
          T1_liquid=290.15,
          Tn_liquid=294.15,
          pstart=176000) annotation (Placement(transformation(
              extent={{-16.5,-14.5},{16.5,14.5}},
              rotation=180,
              origin={17.5,-88.5})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid2(m0=22)
          annotation (Placement(transformation(extent={{60,-116},{80,-96}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid2(T0(
              displayUnit="degC") = 289.15)
          annotation (Placement(transformation(extent={{26,-132},{6,-112}})));
        Components.Controls.sens_PT sens_PT2(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{-4,-94},{-24,-74}})));
        Components.Controls.sens_PT sens_D(redeclare package medium =
              ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
              extent={{6,7},{-6,-7}},
              rotation=-90,
              origin={-67,-60})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT1
          annotation (Placement(transformation(extent={{-102,-8},{-82,12}})));
        Components.Controls.SuperHeat superHeat(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{52,-34},{72,-14}})));
        Components.main_components.Pump_type4 pump_type4_1(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          v_s=1.0395e-5,
          p0=176000,
          pstart_out=1895000,
          T0=295.15,
          v_max=1.9e-5)
          annotation (Placement(transformation(extent={{-58,-106},{-94,-70}})));
        Modelica.Blocks.Sources.Constant const2(k=50)
          annotation (Placement(transformation(extent={{-108,-92},{-98,-82}})));
        Modelica.Blocks.Sources.Constant const1(k=365.65)
          annotation (Placement(transformation(extent={{-112,16},{-102,26}})));
        Components.source_sink_ports.OutFlow_liquid_useM outFlow_liquid_useM
          annotation (Placement(transformation(extent={{-42,14},{-62,34}})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT
          annotation (Placement(transformation(extent={{58,2},{38,22}})));
        Components.Controls.SuperHeat superHeat1(redeclare package medium =
              ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
                extent={{-36,-120},{-56,-100}})));
        Components.main_components.Tank_rou tank_rou(
          L0=0.65,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0=176000,
          T0=296.15,
          Aera=0.5) annotation (Placement(transformation(extent={{-80,-124},{
                  -60,-104}})));
        Modelica.Blocks.Sources.Step step1_5(
          offset=0.2981,
          height=-0.1284,
          startTime=5000)
          annotation (Placement(transformation(extent={{-122,76},{-102,96}})));
        Modelica.Blocks.Sources.Step stepT1_5(
          offset=808.15,
          height=-115,
          startTime=5000)
          annotation (Placement(transformation(extent={{-12,82},{8,102}})));
        Modelica.Blocks.Sources.Step step(height=0.1284, startTime=5600)
          annotation (Placement(transformation(extent={{-122,42},{-102,62}})));
        Modelica.Blocks.Sources.Step stepT1_1(
          height=115,
          offset=0,
          startTime=5600)
          annotation (Placement(transformation(extent={{-14,38},{6,58}})));
        Modelica.Blocks.Math.Add add
          annotation (Placement(transformation(extent={{-80,56},{-60,76}})));
        Modelica.Blocks.Math.Add add1
          annotation (Placement(transformation(extent={{20,62},{40,82}})));
      equation
        connect(Evaporater.OutFlow_working, sens_B.inFlow) annotation (Line(
              points={{13.68,1.52},{21.84,1.52},{21.84,-23.8},{26,-23.8}},
              color={0,0,255}));
        connect(sens_C.outFlow, Evaporater.InFlow_working) annotation (Line(
              points={{-32,-17.8},{-26,-17.8},{-26,1.24},{-18,1.24}}, color={0,
                0,255}));
        connect(preHeater.OutFlow_liquid,outFlow_liquid1. A) annotation (Line(
              points={{-73.84,-44},{-76,-44},{-76,-58},{-84,-58}},
                                                                 color={28,108,
                200}));
        connect(preHeater.OutFlow_working, sens_C.inFlow) annotation (Line(
              points={{-65.48,-24.2},{-65.74,-24.2},{-65.74,-17.8},{-52,-17.8}},
              color={0,0,255}));
        connect(turbine_simple.OutFlow,Condenser. InFlow_working) annotation (
            Line(points={{90.44,-66.22},{90.95,-66.22},{90.95,-83.57},{34,
                -83.57}},
              color={0,0,255}));
        connect(outFlow_liquid2.A,Condenser. OutFlow_liquid) annotation (Line(
              points={{60,-106},{40,-106},{40,-94.88},{34,-94.88}},
                                                                  color={28,108,
                200}));
        connect(sens_PT2.inFlow,Condenser. OutFlow_working) annotation (Line(
              points={{-4,-83.8},{0,-83.8},{0,-83.86},{1.33,-83.86}},color={0,0,
                255}));
        connect(inFlow_liquid2.A,Condenser. InFlow_liquid) annotation (Line(
              points={{6.6,-122.2},{-13.3,-122.2},{-13.3,-95.17},{1.33,-95.17}},
              color={28,108,200}));
        connect(sens_D.outFlow, preHeater.InFlow_working) annotation (Line(
              points={{-67.14,-54},{-65.26,-54},{-65.26,-44}}, color={0,0,255}));
        connect(inFlow_liquid_useT1.A,preHeater. InFlow_liquid) annotation (
            Line(points={{-82.6,1.8},{-82.6,-11.1},{-74.06,-11.1},{-74.06,-24.2}},
              color={28,108,200}));
        connect(sens_B.outFlow, superHeat.inFlow) annotation (Line(points={{46,
                -23.8},{50,-23.8},{52,-23.8}}, color={0,0,255}));
        connect(superHeat.outFlow, turbine_simple.InFlow) annotation (Line(
              points={{72,-23.8},{76,-23.8},{76,-47.86},{77,-47.86}}, color={0,
                0,255}));
        connect(pump_type4_1.OutFlow, sens_D.inFlow) annotation (Line(points={{
                -81.76,-79.72},{-73.88,-79.72},{-73.88,-66},{-67.14,-66}},
              color={0,0,255}));
        connect(const2.y, pump_type4_1.u) annotation (Line(points={{-97.5,-87},
                {-88.75,-87},{-88.75,-85.84},{-82.12,-85.84}}, color={0,0,127}));
        connect(const1.y, inFlow_liquid_useT1.u) annotation (Line(points={{
                -101.5,21},{-101.5,21.5},{-93.2,21.5},{-93.2,6}}, color={0,0,
                127}));
        connect(inFlow_liquid_useT.A, Evaporater.InFlow_liquid) annotation (
            Line(points={{38.6,11.8},{24.3,11.8},{24.3,12.44},{13.68,12.44}},
              color={28,108,200}));
        connect(outFlow_liquid_useM.A, Evaporater.OutFlow_liquid) annotation (
            Line(points={{-42,24},{-32,24},{-32,12.16},{-18,12.16}}, color={28,
                108,200}));
        connect(sens_PT2.outFlow, superHeat1.inFlow) annotation (Line(points={{
                -24,-83.8},{-30,-83.8},{-30,-109.8},{-36,-109.8}}, color={0,0,
                255}));
        connect(tank_rou.InFlow, superHeat1.outFlow) annotation (Line(points={{
                -64,-110},{-56,-110},{-56,-109.8}}, color={0,0,255}));
        connect(pump_type4_1.InFlow, tank_rou.OutFlow) annotation (Line(points=
                {{-68.44,-92.32},{-68.44,-106.16},{-75.4,-106.16},{-75.4,-119}},
              color={0,0,255}));
        connect(step1_5.y, add.u1) annotation (Line(points={{-101,86},{-90.5,86},
                {-90.5,72},{-82,72}}, color={0,0,127}));
        connect(step.y, add.u2) annotation (Line(points={{-101,52},{-92,52},{
                -92,60},{-82,60}}, color={0,0,127}));
        connect(stepT1_5.y, add1.u1) annotation (Line(points={{9,92},{16,92},{
                16,78},{18,78}}, color={0,0,127}));
        connect(stepT1_1.y, add1.u2) annotation (Line(points={{7,48},{16,48},{
                16,66},{18,66}}, color={0,0,127}));
        connect(add.y, outFlow_liquid_useM.u) annotation (Line(points={{-59,66},
                {-52,66},{-52,28.2},{-51.8,28.2}}, color={0,0,127}));
        connect(add1.y, inFlow_liquid_useT.u) annotation (Line(points={{41,72},
                {49.2,72},{49.2,16}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-54,-30},{42,-74}},
                lineColor={28,108,200},
                textString="内燃机余热回收")}));
      end Cycle_condition_change_ORC;

      model Cycle_condition_variation_OSORC
        Components.main_components.Hx1counter_liquid Gas_oil_Hx(
          Cp1=2.34e3,
          U1=1300,
          U2=800,
          Cp_wall=400,
          F_liquid1=3,
          F_liquid2=3,
          V_liquid1=0.7,
          V_liquid2=0.7,
          N=20,
          T1_liquid(displayUnit="degC") = 447.15,
          T1n_liquid(displayUnit="degC") = 499.55,
          T2_liquid(displayUnit="degC") = 808.15,
          T2n_liquid(displayUnit="degC") = 455.42,
          Cp2=1.109e3,
          M2=0.5,
          M_wall=10,
          M1=10)
          annotation (Placement(transformation(extent={{-32,44},{4,70}})));
        Components.main_components.Hx1counter Evaporater(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Utp=3000,
          Uv=800,
          U=1000,
          Cp=2.34e3,
          M_wall=20,
          Cp_wall=400,
          F_working=2.3,
          F_liquid=2.3,
          m_liquid=1,
          N=20,
          m_working=0.7,
          M=10,
          Tin_start=365.33,
          Tout_start=393.15,
          T1_liquid=499.55,
          Tn_liquid=447.85,
          pstart=1895000)
          annotation (Placement(transformation(extent={{-22,-8},{10,20}})));
        Components.main_components.Pump_liquid pump_liquid
          annotation (Placement(transformation(extent={{42,4},{22,24}})));
        Modelica.Blocks.Sources.Constant const(k=1)
          annotation (Placement(transformation(extent={{8,24},{18,34}})));
        Components.Controls.sens_PT sens_B(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{24,-14},{44,6}})));
        Components.Controls.sens_PT sens_C(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{-54,-8},{-34,12}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(m0=2.1)
          annotation (Placement(transformation(extent={{-86,-48},{-106,-28}})));
        Components.main_components.Hx1counter preHeater(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Utp=3000,
          Uv=800,
          U=1000,
          M_wall=20,
          Cp_wall=400,
          N=20,
          M=700,
          Cp=4.2e3,
          m_liquid=2.1,
          F_working=10,
          F_liquid=10,
          m_working=0.7,
          Tin_start=296.45,
          Tout_start=365.35,
          T1_liquid=366.15,
          Tn_liquid=358.45,
          pstart=1895000) annotation (Placement(transformation(
              extent={{-10,-11},{10,11}},
              rotation=90,
              origin={-71,-14})));
        Components.main_components.Turbine_simple turbine_simple(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          Np=35,
          V_s=1.904e-4,
          pin_start=1895000,
          pout_start=176000,
          Tin_start=398.15)
          annotation (Placement(transformation(extent={{68,-52},{96,-18}})));
        Components.main_components.Hx1counter Condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Utp=3000,
          Uv=800,
          U=1000,
          M_wall=20,
          Cp_wall=400,
          N=20,
          Cp=4.2e3,
          m_working=0.7,
          m_liquid=10,
          F_working=20,
          F_liquid=20,
          M=70,
          Tin_start=347.11,
          Tout_start=295.38,
          T1_liquid=290.15,
          Tn_liquid=294.15,
          pstart=176000) annotation (Placement(transformation(
              extent={{-16.5,-14.5},{16.5,14.5}},
              rotation=180,
              origin={15.5,-68.5})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid2(m0=15)
          annotation (Placement(transformation(extent={{58,-96},{78,-76}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid2(T0(
              displayUnit="degC") = 290.15)
          annotation (Placement(transformation(extent={{24,-112},{4,-92}})));
        Components.Controls.sens_PT sens_PT2(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{-6,-74},{-26,-54}})));
        Components.main_components.Tank tank(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          Aera=0.065,
          L0=0.65,
          p0=176000,
          T0=296.15)
          annotation (Placement(transformation(extent={{-60,-84},{-32,-52}})));
        Components.Controls.sens_PT sens_D(redeclare package medium =
              ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
              extent={{6,7},{-6,-7}},
              rotation=-90,
              origin={-69,-40})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT
          annotation (Placement(transformation(extent={{62,54},{42,74}})));
        Components.main_components.Tank_liquid tank_liquid(
          rou=800,
          Cp=2.34e3,
          Area=0.05,
          L0=0.5,
          T0=499.15)
          annotation (Placement(transformation(extent={{52,24},{72,44}})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT1
          annotation (Placement(transformation(extent={{-104,12},{-84,32}})));
        Components.source_sink_ports.OutFlow_liquid_useM outFlow_liquid_useM
          annotation (Placement(transformation(extent={{-72,72},{-52,92}})));
        Components.Controls.SuperHeat superHeat(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{50,-14},{70,6}})));
        Components.main_components.Pump_type4 pump_type4_1(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          v_s=1.0395e-5,
          p0=176000,
          pstart_out=1895000,
          T0=295.15,
          v_max=1.9e-5)
          annotation (Placement(transformation(extent={{-60,-86},{-96,-50}})));
        Modelica.Blocks.Sources.Constant const2(k=50)
          annotation (Placement(transformation(extent={{-110,-72},{-100,-62}})));
        Modelica.Blocks.Sources.Constant const1(k=365.65)
          annotation (Placement(transformation(extent={{-114,36},{-104,46}})));
        Modelica.Blocks.Sources.Step step1_5(
          offset=0.2981,
          height=-0.1284,
          startTime=6000) annotation (Placement(transformation(extent={{-140,
                  136},{-120,156}})));
        Modelica.Blocks.Sources.Step stepT1_5(
          offset=808.15,
          height=-115,
          startTime=6000)
          annotation (Placement(transformation(extent={{0,140},{20,160}})));
        Modelica.Blocks.Sources.Step step(height=0.1284, startTime=6300)
          annotation (Placement(transformation(extent={{-140,102},{-120,122}})));
        Modelica.Blocks.Sources.Step stepT1_1(
          height=115,
          offset=0,
          startTime=6300)
          annotation (Placement(transformation(extent={{-2,96},{18,116}})));
        Modelica.Blocks.Math.Add add
          annotation (Placement(transformation(extent={{-98,114},{-78,134}})));
        Modelica.Blocks.Math.Add add1
          annotation (Placement(transformation(extent={{32,120},{52,140}})));
        Modelica.Blocks.Sources.Constant const3(k=0.2981)
          annotation (Placement(transformation(extent={{-100,68},{-90,78}})));
        Modelica.Blocks.Sources.Constant const4(k=808.15)
          annotation (Placement(transformation(extent={{78,64},{88,74}})));
        Modelica.Blocks.Sources.Pulse pulse(
          offset=0.2981,
          amplitude=-0.1284,
          period=60,
          nperiod=3,
          startTime=2000)
          annotation (Placement(transformation(extent={{-66,174},{-46,194}})));
      equation
        connect(Gas_oil_Hx.inflow1,Evaporater. OutFlow_liquid) annotation (Line(
              points={{-32.36,51.8},{-32.36,50.9},{-44,50.9},{-44,12.16},{-22,
                12.16}}, color={28,108,200}));
        connect(Evaporater.InFlow_liquid,pump_liquid. outFlow) annotation (Line(
              points={{9.68,12.44},{16.84,12.44},{16.84,14},{22,14}},  color={
                28,108,200}));
        connect(const.y,pump_liquid. u) annotation (Line(points={{18.5,29},{
                18.5,28.5},{31.9,28.5},{31.9,17.9}}, color={0,0,127}));
        connect(Evaporater.OutFlow_working, sens_B.inFlow) annotation (Line(
              points={{9.68,1.52},{19.84,1.52},{19.84,-3.8},{24,-3.8}}, color={
                0,0,255}));
        connect(sens_C.outFlow, Evaporater.InFlow_working) annotation (Line(
              points={{-34,2.2},{-28,2.2},{-28,1.24},{-22,1.24}}, color={0,0,
                255}));
        connect(preHeater.OutFlow_liquid,outFlow_liquid1. A) annotation (Line(
              points={{-75.84,-24},{-78,-24},{-78,-38},{-86,-38}},
                                                                 color={28,108,
                200}));
        connect(preHeater.OutFlow_working, sens_C.inFlow) annotation (Line(
              points={{-67.48,-4.2},{-67.74,-4.2},{-67.74,2.2},{-54,2.2}},
              color={0,0,255}));
        connect(turbine_simple.OutFlow,Condenser. InFlow_working) annotation (
            Line(points={{88.44,-46.22},{88.95,-46.22},{88.95,-63.57},{32,
                -63.57}},
              color={0,0,255}));
        connect(outFlow_liquid2.A,Condenser. OutFlow_liquid) annotation (Line(
              points={{58,-86},{38,-86},{38,-74.88},{32,-74.88}}, color={28,108,
                200}));
        connect(sens_PT2.inFlow,Condenser. OutFlow_working) annotation (Line(
              points={{-6,-63.8},{-2,-63.8},{-2,-63.86},{-0.67,-63.86}},
                                                                     color={0,0,
                255}));
        connect(inFlow_liquid2.A,Condenser. InFlow_liquid) annotation (Line(
              points={{4.6,-102.2},{-15.3,-102.2},{-15.3,-75.17},{-0.67,-75.17}},
              color={28,108,200}));
        connect(tank.InFlow,sens_PT2. outFlow) annotation (Line(points={{-37.6,
                -61.6},{-38,-61.6},{-38,-63.8},{-26,-63.8}}, color={0,0,255}));
        connect(sens_D.outFlow, preHeater.InFlow_working) annotation (Line(
              points={{-69.14,-34},{-67.26,-34},{-67.26,-24}}, color={0,0,255}));
        connect(inFlow_liquid_useT.A,Gas_oil_Hx. inflow2) annotation (Line(
              points={{42.6,63.8},{18.3,63.8},{18.3,63.24},{3.64,63.24}},
              color={28,108,200}));
        connect(Gas_oil_Hx.outflow1,tank_liquid. inFlow) annotation (Line(
              points={{3.64,51.8},{72.82,51.8},{72.82,38.2},{68.2,38.2}},
              color={28,108,200}));
        connect(tank_liquid.outFlow,pump_liquid. inFlow) annotation (Line(
              points={{57.2,28.4},{57.2,13.2},{42,13.2},{42,14}}, color={28,108,
                200}));
        connect(inFlow_liquid_useT1.A,preHeater. InFlow_liquid) annotation (
            Line(points={{-84.6,21.8},{-84.6,8.9},{-76.06,8.9},{-76.06,-4.2}},
              color={28,108,200}));
        connect(outFlow_liquid_useM.A,Gas_oil_Hx. outflow2) annotation (Line(
              points={{-72,82},{-76,82},{-76,63.5},{-32,63.5}},   color={28,108,
                200}));
        connect(sens_B.outFlow, superHeat.inFlow) annotation (Line(points={{44,
                -3.8},{48,-3.8},{50,-3.8}}, color={0,0,255}));
        connect(superHeat.outFlow, turbine_simple.InFlow) annotation (Line(
              points={{70,-3.8},{74,-3.8},{74,-27.86},{75,-27.86}}, color={0,0,
                255}));
        connect(pump_type4_1.InFlow, tank.OutFlow) annotation (Line(points={{
                -70.44,-72.32},{-61.22,-72.32},{-61.22,-76},{-53.56,-76}},
              color={0,0,255}));
        connect(pump_type4_1.OutFlow, sens_D.inFlow) annotation (Line(points={{
                -83.76,-59.72},{-75.88,-59.72},{-75.88,-46},{-69.14,-46}},
              color={0,0,255}));
        connect(const2.y, pump_type4_1.u) annotation (Line(points={{-99.5,-67},
                {-90.75,-67},{-90.75,-65.84},{-84.12,-65.84}}, color={0,0,127}));
        connect(const1.y, inFlow_liquid_useT1.u) annotation (Line(points={{
                -103.5,41},{-103.5,41.5},{-95.2,41.5},{-95.2,26}}, color={0,0,
                127}));
        connect(step1_5.y, add.u1) annotation (Line(points={{-119,146},{-108.5,
                146},{-108.5,130},{-100,130}}, color={0,0,127}));
        connect(step.y, add.u2) annotation (Line(points={{-119,112},{-110,112},
                {-110,118},{-100,118}}, color={0,0,127}));
        connect(stepT1_5.y, add1.u1) annotation (Line(points={{21,150},{28,150},
                {28,136},{30,136}}, color={0,0,127}));
        connect(stepT1_1.y, add1.u2) annotation (Line(points={{19,106},{28,106},
                {28,124},{30,124}}, color={0,0,127}));
        connect(add.y, outFlow_liquid_useM.u) annotation (Line(points={{-77,124},
                {-70,124},{-70,86.2},{-62.2,86.2}}, color={0,0,127}));
        connect(add1.y, inFlow_liquid_useT.u) annotation (Line(points={{53,130},
                {53.2,130},{53.2,68}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-56,-10},{40,-54}},
                lineColor={28,108,200},
                textString="内燃机余热回收")}));
      end Cycle_condition_variation_OSORC;
    end TestCycle2;

    package TestMedium

      model TestUSDr125
        parameter Modelica.SIunits.Temperature T=320;
        parameter Modelica.SIunits.Pressure p0=1.4463e6;
        parameter Modelica.SIunits.SpecificEnthalpy h0=2.2e5;
       TJUthermo.DataRecord.SatProperty Sat1;
        Modelica.SIunits.SpecificEnthalpy h;
        TJUthermo.Media.MediaUSD.R125_USD.R125state_ph state(P=p0,h=h0);
      equation
        Sat1=TJUthermo.Media.MediaUSD.R125_USD.Sat_property_T(T);
        h=TJUthermo.Media.MediaUSD.R125_USD.SpecificEnthalpy_Trou(T,500);
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestUSDr125;

      model TestR421mix
      extends Modelica.Icons.Example;
        replaceable package Medium=TJUthermo.Media.R421mix;
        Real X[Medium.nC]=Medium.reference_X;
        Real MMX[Medium.nC]=Medium.MMX;
        Real MM;
       //Modelica.SIunits.Temperature T;
       Medium.Auxilary.Properties properties;
      //MultiPhaseMixture.Interfaces.FlashProperties  equ;
        Modelica.SIunits.Pressure P;//(start=0.8e6);
        parameter Modelica.SIunits.SpecificEnthalpy h=3.5e5;
        parameter Modelica.SIunits.MassFlowRate m1=1;
        parameter Modelica.SIunits.MassFlowRate m2=1.08;
        parameter Modelica.SIunits.Volume V=0.004;
        Modelica.SIunits.Density rou;//(start=50);
        Real roudh;
        Real roudP(start=0.7e-4);
        Real dt;
        Real roud;
        Real tt;
      equation
       MM= Medium.Auxilary.getAverageMolarMass(X_unit=1,eo=Medium.eo);
       // State=MultiPhaseMixture.PreDefined.Mixtures.RefProp.R421A.Auxilary.calcProperties_phX(P,h,X_unit=1,eo=Medium.eo);
       properties=Medium.Auxilary.calcProperties_phX(p=P,h=h,X=X,X_unit=1,phase=0,eo=Medium.eo);
       //T=Medium.Auxilary.calcProperties_phX(P,h,X,eo=Medium.eo);
       //equ=Medium.Auxilary.Wrapper.equilibrium_X(X,properties);
       der(rou)=roudP*der(P);
       //2*der(P)/(properties.dpdd_TN_1ph[1]+properties.dpdd_TN_1ph[2]);
       m2-m1=V*der(rou);
        roudh=Medium.Auxilary.Wrapper.density_derh_p_X(X=X,state=properties,eo=Medium.eo);
        //if time
        dt=floor(time);
        roud=Medium.Auxilary.Wrapper.density_derp_h_X(X=X,state=properties,eo=Medium.eo);
        tt=time;
        when abs(time-dt)<0.01 then
          roudP=roud;
        end when;
      initial equation
        P=0.8e6;
        rou=50;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestR421mix;

      model TestTolune
        replaceable package medium = ExternalMedia.Media.CoolPropMedium(mediumName="toluene");
          parameter Modelica.SIunits.Pressure p0=1.1e5;
          Modelica.SIunits.Temperature T;
          parameter Modelica.SIunits.SpecificEnthalpy h=-49338;//-114608;//-4.9e-5;//-490000;//-496538;
          medium.ThermodynamicState state;
          medium.SaturationProperties Sat1;
      equation
           state=medium.setState_ph(p0,h);
           Sat1=medium.setSat_p(p0);
           T=state.T;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestTolune;

      model TestTolune_2
      replaceable package medium = ExternalMedia.Media.CoolPropMedium(mediumName="toluene");
          parameter Modelica.SIunits.Pressure p0=1.1e5;
          parameter Modelica.SIunits.Temperature T=273.15+50;
          Modelica.SIunits.SpecificEnthalpy h;//=-22616;//-4.9e-5;//-490000;//-496538;
          medium.ThermodynamicState state;
          medium.SaturationProperties Sat1;
      equation
           state=medium.setState_pT(p0,T);
           Sat1=medium.setSat_p(p0);
           h=state.h;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestTolune_2;
    end TestMedium;

    package Validationsystem_taiwan

      model validation1
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          N=20,
          V_working=0.1753 + 1.1414,
          V_liquid=0.1224 + 0.2279,
          M=0.4*1000,
          Cp=4.2e3,
          Cp_wall=390,
          m_working=11.85,
          m_liquid=26.7,
          pstart(displayUnit="kPa") = 1440000,
          M_wall=50,
          F_working=123,
          F_liquid=123,
          Ul=2350,
          Utp=11000,
          Uv=1700,
          U=3500,
          Tin_start=315.15,
          Tout_start=383.15,
          T1_liquid=393.15,
          Tn_liquid=368.15)
          annotation (Placement(transformation(extent={{-44,14},{4,64}})));
        Components.main_components.Hx1counter hx1counter1(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          N=20,
          V_working=0.4925,
          V_liquid=0.9611,
          M=0.9611*1000,
          Cp=4200,
          M_wall=55,
          Cp_wall=390,
          m_working=11.85,
          m_liquid=30,
          pstart(displayUnit="kPa") = 250000,
          F_working=150,
          F_liquid=150,
          Ul=2300,
          Utp=8500,
          Uv=1300,
          U=3000,
          Tin_start=316.15,
          Tout_start=293.15,
          T1_liquid=288.15,
          Tn_liquid=299.15)
          annotation (Placement(transformation(extent={{42,-36},{-6,-80}})));
        Components.main_components.Tank tank(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0(displayUnit="kPa") = 250460,
          Aera=1.2,
          L0=1,
          T0=303.15)
          annotation (Placement(transformation(extent={{-58,-54},{-12,-8}})));
        Components.main_components.Turbine_simple turbine_simple(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          n_turbine=0.681,
          Np=55,
          V_s=0.0025,
          n_ele=0.89,
          pin_start(displayUnit="kPa") = 1440000,
          pout_start(displayUnit="kPa") = 250000,
          Tin_start(displayUnit="degC") = 383.15)
          annotation (Placement(transformation(extent={{16,-12},{60,32}})));
        Components.main_components.pump_type3 pump_type3_1(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          n_pump=0.681,
          p0(displayUnit="kPa") = 249226,
          pstart_out(displayUnit="kPa") = 1440000,
          T0=310.15,
          v_s=1.7423e-4)
          annotation (Placement(transformation(extent={{-38,-30},{-78,20}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid(m0=26.7)
          annotation (Placement(transformation(extent={{-58,70},{-38,90}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid1(T0=298.15)
          annotation (Placement(transformation(extent={{-62,-84},{-38,-62}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(m0=30)
          annotation (Placement(transformation(extent={{60,-80},{80,-60}})));
        Modelica.Blocks.Sources.Constant const(k=50)
          annotation (Placement(transformation(extent={{-98,-14},{-78,6}})));
        Components.Controls.sens_PT sens_PT(redeclare package medium =
              ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-68,-34})));
        Components.Controls.sens_PT sens_PT1(redeclare package medium =
              ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=90,
              origin={56,-22})));
        Components.Controls.sens_PT sens_PT2(redeclare package medium =
              ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-6,-26})));
        Components.Controls.SuperHeat superHeat(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{12,32},{26,48}})));
        Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(
          tableOnFile=true,
          tableName="source1",
          fileName="validation_s.mat",
          columns={2})
          annotation (Placement(transformation(extent={{22,82},{34,94}})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT
          annotation (Placement(transformation(extent={{94,42},{74,62}})));
        Modelica.Blocks.Math.Add add
          annotation (Placement(transformation(extent={{56,68},{76,88}})));
        Modelica.Blocks.Sources.Constant const1(k=273.15)
          annotation (Placement(transformation(extent={{22,64},{34,76}})));
      equation
        connect(const.y, pump_type3_1.u) annotation (Line(points={{-77,-4},{-70,
                -4},{-70,-2},{-64,-2}}, color={0,0,127}));
        connect(outFlow_liquid.A, hx1counter.OutFlow_liquid) annotation (Line(
              points={{-58,80},{-70,80},{-70,50},{-44,50}}, color={28,108,200}));
        connect(inFlow_liquid1.A, hx1counter1.InFlow_liquid) annotation (Line(
              points={{-38.72,-73.22},{-25.36,-73.22},{-25.36,-68.12},{-5.52,
                -68.12}}, color={28,108,200}));
        connect(hx1counter1.OutFlow_liquid, outFlow_liquid1.A) annotation (Line(
              points={{42,-67.68},{50,-67.68},{50,-70},{60,-70}}, color={28,108,
                200}));
        connect(pump_type3_1.OutFlow, hx1counter.InFlow_working) annotation (
            Line(points={{-64.4,6.5},{-64.4,18.25},{-44,18.25},{-44,30.5}},
              color={0,0,255}));
        connect(tank.OutFlow, sens_PT.inFlow) annotation (Line(points={{-47.42,
                -42.5},{-56.71,-42.5},{-56.71,-44},{-68.2,-44}}, color={0,0,255}));
        connect(sens_PT.outFlow, pump_type3_1.InFlow) annotation (Line(points={
                {-68.2,-24},{-50,-24},{-50,-11},{-49.6,-11}}, color={0,0,255}));
        connect(turbine_simple.OutFlow, sens_PT1.inFlow) annotation (Line(
              points={{48.12,-4.52},{48.12,-11.26},{55.8,-11.26},{55.8,-12}},
              color={0,0,255}));
        connect(sens_PT1.outFlow, hx1counter1.InFlow_working) annotation (Line(
              points={{55.8,-32},{46,-32},{46,-50.52},{42,-50.52}}, color={0,0,
                255}));
        connect(tank.InFlow, sens_PT2.outFlow) annotation (Line(points={{-21.2,
                -21.8},{-14.6,-21.8},{-14.6,-16},{-6.2,-16}}, color={0,0,255}));
        connect(hx1counter1.OutFlow_working, sens_PT2.inFlow) annotation (Line(
              points={{-5.52,-50.96},{-5.52,-43.48},{-6.2,-43.48},{-6.2,-36}},
              color={0,0,255}));
        connect(hx1counter.OutFlow_working, superHeat.inFlow) annotation (Line(
              points={{3.52,31},{10.76,31},{10.76,40.16},{12,40.16}}, color={0,
                0,255}));
        connect(superHeat.outFlow, turbine_simple.InFlow) annotation (Line(
              points={{26,40.16},{27,40.16},{27,19.24}}, color={0,0,255}));
        connect(add.y, inFlow_liquid_useT.u) annotation (Line(points={{77,78},{
                86,78},{86,56},{85.2,56}}, color={0,0,127}));
        connect(combiTimeTable.y[1], add.u1) annotation (Line(points={{34.6,88},
                {46,88},{46,84},{54,84}}, color={0,0,127}));
        connect(const1.y, add.u2) annotation (Line(points={{34.6,70},{42,70},{
                42,72},{54,72}}, color={0,0,127}));
        connect(inFlow_liquid_useT.A, hx1counter.InFlow_liquid) annotation (
            Line(points={{74.6,51.8},{38.3,51.8},{38.3,50.5},{3.52,50.5}},
              color={28,108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end validation1;

      model validation_with_pump3
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          N=20,
          V_working=0.1753 + 1.1414,
          V_liquid=0.1224 + 0.2279,
          M=0.4*1000,
          Cp=4.2e3,
          Cp_wall=390,
          m_working=11.85,
          m_liquid=26.7,
          pstart(displayUnit="kPa") = 1440000,
          M_wall=50,
          F_working=123,
          F_liquid=123,
          Ul=2350,
          Utp=11000,
          Uv=1700,
          U=3500,
          Tin_start=315.15,
          Tout_start=383.15,
          T1_liquid=393.15,
          Tn_liquid=368.15)
          annotation (Placement(transformation(extent={{-34,24},{14,74}})));
        Components.main_components.Hx1counter hx1counter1(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          N=20,
          V_working=0.4925,
          V_liquid=0.9611,
          M=0.9611*1000,
          Cp=4200,
          M_wall=55,
          Cp_wall=390,
          m_working=11.85,
          m_liquid=30,
          pstart(displayUnit="kPa") = 250000,
          F_working=150,
          F_liquid=150,
          Ul=2300,
          Utp=8500,
          Uv=1300,
          U=3000,
          Tin_start=316.15,
          Tout_start=293.15,
          T1_liquid=288.15,
          Tn_liquid=299.15)
          annotation (Placement(transformation(extent={{52,-26},{4,-70}})));
        Components.main_components.Tank tank(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0(displayUnit="kPa") = 250460,
          Aera=1.2,
          L0=1,
          T0=303.15)
          annotation (Placement(transformation(extent={{-48,-44},{-2,2}})));
        Components.main_components.Turbine_simple turbine_simple(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          n_turbine=0.681,
          Np=55,
          V_s=0.0025,
          n_ele=0.89,
          pin_start(displayUnit="kPa") = 1440000,
          pout_start(displayUnit="kPa") = 250000,
          Tin_start(displayUnit="degC") = 383.15)
          annotation (Placement(transformation(extent={{26,-2},{70,42}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid(m0=26.7)
          annotation (Placement(transformation(extent={{-48,80},{-28,100}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid1(T0=298.15)
          annotation (Placement(transformation(extent={{-52,-74},{-28,-52}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(m0=30)
          annotation (Placement(transformation(extent={{70,-70},{90,-50}})));
        Modelica.Blocks.Sources.Constant const(k=50)
          annotation (Placement(transformation(extent={{-88,-4},{-68,16}})));
        Components.Controls.sens_PT sens_PT(redeclare package medium =
              ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-58,-24})));
        Components.Controls.sens_PT sens_PT1(redeclare package medium =
              ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=90,
              origin={66,-12})));
        Components.Controls.sens_PT sens_PT2(redeclare package medium =
              ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={4,-16})));
        Components.Controls.SuperHeat superHeat(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{22,42},{36,58}})));
        Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(
          tableOnFile=true,
          tableName="source1",
          fileName="validation_s.mat",
          columns={2})
          annotation (Placement(transformation(extent={{32,92},{44,104}})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT
          annotation (Placement(transformation(extent={{104,52},{84,72}})));
        Modelica.Blocks.Math.Add add
          annotation (Placement(transformation(extent={{66,78},{86,98}})));
        Modelica.Blocks.Sources.Constant const1(k=273.15)
          annotation (Placement(transformation(extent={{32,74},{44,86}})));
        Components.main_components.Pump_type4 pump_type4_1(
          v_s=1.7423e-4,
          p0(displayUnit="kPa") = 249226,
          pstart_out(displayUnit="kPa") = 1440000,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          T0=310.15,
          v_max=0.0001)
          annotation (Placement(transformation(extent={{-28,-4},{-64,32}})));
      equation
        connect(outFlow_liquid.A,hx1counter. OutFlow_liquid) annotation (Line(
              points={{-48,90},{-60,90},{-60,60},{-34,60}}, color={28,108,200}));
        connect(inFlow_liquid1.A,hx1counter1. InFlow_liquid) annotation (Line(
              points={{-28.72,-63.22},{-15.36,-63.22},{-15.36,-58.12},{4.48,
                -58.12}}, color={28,108,200}));
        connect(hx1counter1.OutFlow_liquid,outFlow_liquid1. A) annotation (Line(
              points={{52,-57.68},{60,-57.68},{60,-60},{70,-60}}, color={28,108,
                200}));
        connect(tank.OutFlow,sens_PT. inFlow) annotation (Line(points={{-37.42,
                -32.5},{-46.71,-32.5},{-46.71,-34},{-58.2,-34}}, color={0,0,255}));
        connect(turbine_simple.OutFlow,sens_PT1. inFlow) annotation (Line(
              points={{58.12,5.48},{58.12,-1.26},{65.8,-1.26},{65.8,-2}},
              color={0,0,255}));
        connect(sens_PT1.outFlow,hx1counter1. InFlow_working) annotation (Line(
              points={{65.8,-22},{56,-22},{56,-40.52},{52,-40.52}}, color={0,0,
                255}));
        connect(tank.InFlow,sens_PT2. outFlow) annotation (Line(points={{-11.2,
                -11.8},{-4.6,-11.8},{-4.6,-6},{3.8,-6}},      color={0,0,255}));
        connect(hx1counter1.OutFlow_working,sens_PT2. inFlow) annotation (Line(
              points={{4.48,-40.96},{4.48,-33.48},{3.8,-33.48},{3.8,-26}},
              color={0,0,255}));
        connect(hx1counter.OutFlow_working,superHeat. inFlow) annotation (Line(
              points={{13.52,41},{20.76,41},{20.76,50.16},{22,50.16}},color={0,
                0,255}));
        connect(superHeat.outFlow,turbine_simple. InFlow) annotation (Line(
              points={{36,50.16},{37,50.16},{37,29.24}}, color={0,0,255}));
        connect(add.y,inFlow_liquid_useT. u) annotation (Line(points={{87,88},{
                96,88},{96,66},{95.2,66}}, color={0,0,127}));
        connect(combiTimeTable.y[1],add. u1) annotation (Line(points={{44.6,98},
                {56,98},{56,94},{64,94}}, color={0,0,127}));
        connect(const1.y,add. u2) annotation (Line(points={{44.6,80},{52,80},{
                52,82},{64,82}}, color={0,0,127}));
        connect(inFlow_liquid_useT.A,hx1counter. InFlow_liquid) annotation (
            Line(points={{84.6,61.8},{48.3,61.8},{48.3,60.5},{13.52,60.5}},
              color={28,108,200}));
        connect(const.y, pump_type4_1.u) annotation (Line(points={{-67,6},{
                -61.5,6},{-61.5,16.16},{-52.12,16.16}}, color={0,0,127}));
        connect(pump_type4_1.OutFlow, hx1counter.InFlow_working) annotation (
            Line(points={{-51.76,22.28},{-42.88,22.28},{-42.88,40.5},{-34,40.5}},
              color={0,0,255}));
        connect(pump_type4_1.InFlow, sens_PT.outFlow) annotation (Line(points={
                {-38.44,9.68},{-38.44,-2.16},{-58.2,-2.16},{-58.2,-14}}, color=
                {0,0,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end validation_with_pump3;
    end Validationsystem_taiwan;

    package TestToluene_Cycle
      model SolarCycle_1
        Components.main_components.Hx1counter hx1counter1(
          V_working=0.4925,
          V_liquid=0.9611,
          M=0.9611*1000,
          Cp=4200,
          M_wall=55,
          Cp_wall=390,
          redeclare package Medium = ThermoCycle.Media.Toluene_CP,
          Tin_start(displayUnit="K") = 424,
          m_working=5,
          m_liquid=25,
          T1_liquid(displayUnit="K") = 348.15,
          Tn_liquid(displayUnit="K") = 368.15,
          witdth_x=0.05,
          Ul=2100,
          Utp=2300,
          Uv=1500,
          U=1100,
          Tout_start(displayUnit="K") = 360,
          F_working=100,
          F_liquid=100,
          pstart(displayUnit="kPa") = 110000,
          N=10,
          max_drhode=50)
          annotation (Placement(transformation(extent={{20,30},{-28,-14}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid1(T0(
              displayUnit="K") = 370)
          annotation (Placement(transformation(extent={{-76,-24},{-52,-2}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(m0=10)
          annotation (Placement(transformation(extent={{48,-26},{28,-6}})));
        Components.source_sink_ports.OutFlow outFlow2(
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          m0=5,
          p0(displayUnit="kPa") = 102550)
          annotation (Placement(transformation(extent={{-84,18},{-64,38}})));
        Components.source_sink_ports.InFlow inFlow2(
          UseT=true,
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          m0=-5,
          p0(displayUnit="kPa") = 102550,
          T(displayUnit="K") = 424.3)
          annotation (Placement(transformation(extent={{68,20},{88,40}})));
        Components.Controls.sens_PT sens_PT(redeclare package medium =
              ThermoCycle.Media.Toluene_CP)
          annotation (Placement(transformation(extent={{-46,2},{-66,22}})));
        Components.Controls.sens_PT sens_PT1(
                                            redeclare package medium =
              ThermoCycle.Media.Toluene_CP)
          annotation (Placement(transformation(extent={{58,4},{38,24}})));
      equation
        connect(inFlow_liquid1.A,hx1counter1. InFlow_liquid) annotation (Line(
              points={{-52.72,-13.22},{-39.36,-13.22},{-39.36,-2.12},{-27.52,
                -2.12}},  color={28,108,200}));
        connect(hx1counter1.OutFlow_liquid,outFlow_liquid1. A) annotation (Line(
              points={{20,-1.68},{54,-1.68},{54,-16},{48,-16}},   color={28,108,
                200}));
        connect(outFlow2.node, sens_PT.outFlow) annotation (Line(points={{-65.8,
                27.8},{-65.8,22.9},{-66,22.9},{-66,12.2}}, color={0,0,255}));
        connect(sens_PT.inFlow, hx1counter1.OutFlow_working) annotation (Line(
              points={{-46,12.2},{-38,12.2},{-38,15.04},{-27.52,15.04}}, color=
                {0,0,255}));
        connect(sens_PT1.inFlow, inFlow2.node) annotation (Line(points={{58,14.2},
                {65,14.2},{65,29.6},{70.4,29.6}},       color={0,0,255}));
        connect(sens_PT1.outFlow, hx1counter1.InFlow_working) annotation (Line(
              points={{38,14.2},{26,14.2},{26,15.48},{20,15.48}}, color={0,0,
                255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end SolarCycle_1;

    public
      model testModel
        Components.main_components.pump_type3 pump_type3_1(
          n_pump=0.681,
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          p0(displayUnit="kPa") = 102500,
          pstart_out(displayUnit="kPa") = 696140,
          T0(displayUnit="K") = 383,
          v_s=1.2867e-4)
          annotation (Placement(transformation(extent={{-26,22},{-66,72}})));
        Modelica.Blocks.Sources.Constant const(k=50)
          annotation (Placement(transformation(extent={{-82,44},{-62,64}})));
        Components.source_sink_ports.OutFlow outFlow(
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          m0=5,
          p0(displayUnit="kPa") = 696140)
          annotation (Placement(transformation(extent={{-84,74},{-64,94}})));
        Components.source_sink_ports.InFlow inFlow1(
          UseT=true,
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          m0=-5,
          p0(displayUnit="kPa") = 102500,
          T(displayUnit="K") = 383)
          annotation (Placement(transformation(extent={{-54,10},{-34,30}})));
        Components.main_components.Hx1counter hx1counter1(
          N=20,
          V_working=0.4925,
          V_liquid=0.9611,
          M=0.9611*1000,
          Cp=4200,
          M_wall=55,
          Cp_wall=390,
          m_working=11.85,
          m_liquid=30,
          Ul=2300,
          Utp=8500,
          Uv=1300,
          U=3000,
          redeclare package Medium = ThermoCycle.Media.Toluene_CP,
          F_working=105.89,
          F_liquid=120.2,
          Tin_start(displayUnit="K") = 424,
          Tout_start(displayUnit="K") = 384,
          T1_liquid=288.15,
          Tn_liquid=299.15,
          pstart(displayUnit="kPa") = 102550)
          annotation (Placement(transformation(extent={{62,-16},{14,-60}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid1(T0=298.15)
          annotation (Placement(transformation(extent={{-20,-64},{4,-42}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(m0=45)
          annotation (Placement(transformation(extent={{80,-60},{100,-40}})));
        Components.source_sink_ports.OutFlow outFlow1(
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          m0=5,
          p0(displayUnit="kPa") = 102500)
          annotation (Placement(transformation(extent={{-20,-42},{0,-22}})));
        Components.main_components.Turbine_simple turbine_simple(
          n_turbine=0.681,
          n_ele=0.89,
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          pin_start(displayUnit="kPa") = 696000,
          pout_start(displayUnit="kPa") = 102000,
          Tin_start(displayUnit="K") = 468.99,
          V_s=0.0051,
          Np=50)
          annotation (Placement(transformation(extent={{28,-2},{72,42}})));
        Components.main_components.Turbine_simple turbine_simple1(
          n_turbine=0.681,
          n_ele=0.89,
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          pin_start(displayUnit="kPa") = 696000,
          pout_start(displayUnit="kPa") = 102000,
          Tin_start(displayUnit="K") = 468.99,
          Np=50,
          V_s=0.0063)
          annotation (Placement(transformation(extent={{-2,36},{42,80}})));
        Components.source_sink_ports.InFlow inFlow2(
          UseT=true,
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          p0(displayUnit="kPa") = 696000,
          m0=-5,
          T(displayUnit="K") = 469)
          annotation (Placement(transformation(extent={{-16,58},{-36,78}})));
        Components.source_sink_ports.OutFlow outFlow2(
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          m0=5,
          p0(displayUnit="kPa") = 102500)
          annotation (Placement(transformation(extent={{-6,22},{14,42}})));
        Components.source_sink_ports.InFlow inFlow3(
          UseT=true,
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          p0(displayUnit="kPa") = 696000,
          m0=-5,
          T(displayUnit="K") = 424)
          annotation (Placement(transformation(extent={{40,-14},{20,6}})));
        Components.main_components.Pump_type2 pump_type2_1
          annotation (Placement(transformation(extent={{-162,-70},{-86,4}})));
        Components.main_components.Turbine_simple turbine_simple2
          annotation (Placement(transformation(extent={{-146,-10},{-92,34}})));
        Components.main_components.Tank tank
          annotation (Placement(transformation(extent={{-108,-92},{-48,-22}})));
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.Toluene_CP,
          M=700,
          Cp=2400,
          M_wall=50,
          Cp_wall=380,
          m_working=5,
          m_liquid=20,
          V_working=1,
          V_liquid=1,
          U=1100,
          N=20,
          max_drhode=10,
          witdth_x=1,
          F_working=160,
          F_liquid=160,
          Tin_start=413.15,
          Tout_start=373.15,
          T1_liquid=343.15,
          Tn_liquid=403.15,
          pstart=110000)
          annotation (Placement(transformation(extent={{154,-74},{78,-128}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid(m0=20)
          annotation (Placement(transformation(extent={{84,-76},{104,-56}})));
        Components.main_components.Tank tank1(
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          p0=110000,
          T0=380.15,
          Aera=3,
          L0=1)
          annotation (Placement(transformation(extent={{10,-112},{42,-80}})));
        Components.source_sink_ports.InFlow inFlow4(
          UseT=true,
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          m0=-5,
          p0(displayUnit="bar") = 110000,
          T(displayUnit="degC") = 393.15)
          annotation (Placement(transformation(extent={{144,-52},{124,-32}})));
      equation
        connect(const.y,pump_type3_1. u) annotation (Line(points={{-61,54},{-54,
                54},{-54,50},{-52,50}}, color={0,0,127}));
        connect(outFlow.node, pump_type3_1.OutFlow) annotation (Line(points={{-65.8,
                83.8},{-53.9,83.8},{-53.9,58.5},{-52.4,58.5}},       color={0,0,
                255}));
        connect(inFlow1.node, pump_type3_1.InFlow) annotation (Line(points={{-51.6,
                19.6},{-51.6,32.8},{-37.6,32.8},{-37.6,41}},       color={0,0,
                255}));
        connect(inFlow_liquid1.A,hx1counter1. InFlow_liquid) annotation (Line(
              points={{3.28,-53.22},{8.64,-53.22},{8.64,-48.12},{14.48,-48.12}},
                          color={28,108,200}));
        connect(hx1counter1.OutFlow_liquid,outFlow_liquid1. A) annotation (Line(
              points={{62,-47.68},{70,-47.68},{70,-50},{80,-50}}, color={28,108,
                200}));
        connect(outFlow1.node, hx1counter1.OutFlow_working) annotation (Line(
              points={{-1.8,-32.2},{2.1,-32.2},{2.1,-30.96},{14.48,-30.96}},
              color={0,0,255}));
        connect(turbine_simple.OutFlow, hx1counter1.InFlow_working) annotation (
           Line(points={{60.12,5.48},{60.12,-7.26},{62,-7.26},{62,-30.52}},
              color={0,0,255}));
        connect(inFlow2.node, turbine_simple1.InFlow) annotation (Line(points={{-18.4,
                67.6},{-12.2,67.6},{-12.2,67.24},{9,67.24}},        color={0,0,
                255}));
        connect(turbine_simple1.OutFlow, outFlow2.node) annotation (Line(points=
               {{30.12,43.48},{22.06,43.48},{22.06,31.8},{12.2,31.8}}, color={0,
                0,255}));
        connect(inFlow_liquid.A, hx1counter.InFlow_liquid) annotation (Line(
              points={{-60.6,-80.2},{-60.6,-74.1},{78.76,-74.1},{78.76,-113.42}},
              color={28,108,200}));
        connect(outFlow_liquid.A, hx1counter.OutFlow_liquid) annotation (Line(
              points={{84,-66},{62,-66},{62,-112.88},{154,-112.88}}, color={28,
                108,200}));
        connect(tank1.InFlow, hx1counter.OutFlow_working) annotation (Line(
              points={{35.6,-89.6},{-44,-89.6},{-44,-92.36},{78.76,-92.36}},
              color={0,0,255}));
        connect(inFlow4.node, hx1counter.InFlow_working) annotation (Line(
              points={{141.6,-42.4},{141.6,-30.2},{154,-30.2},{154,-91.82}},
              color={0,0,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end testModel;

      model testModel2
        Components.main_components.Hx1counter hx1counter(
          N=20,
          V_working=0.1753 + 1.1414,
          V_liquid=0.1224 + 0.2279,
          M=0.4*1000,
          Cp=4.2e3,
          Cp_wall=390,
          M_wall=50,
          Ul=2350,
          Uv=1700,
          redeclare package Medium = ThermoCycle.Media.Toluene_CP,
          F_working=100.05,
          F_liquid=120,
          m_working=5,
          m_liquid=25,
          pstart(displayUnit="kPa") = 696000,
          T1_liquid(displayUnit="K") = 527.91,
          Tn_liquid(displayUnit="K") = 463.65,
          Utp=8000,
          U=3000,
          Tin_start(displayUnit="K") = 383,
          Tout_start(displayUnit="K") = 470)
          annotation (Placement(transformation(extent={{-48,-20},{0,30}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid(T0(displayUnit=
               "K") = 527)
          annotation (Placement(transformation(extent={{34,30},{54,50}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid(m0=15.91)
          annotation (Placement(transformation(extent={{-62,36},{-42,56}})));
        Components.source_sink_ports.InFlow inFlow(
          UseT=true,
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          p0(displayUnit="kPa") = 696000,
          T(displayUnit="K") = 389,
          m0=-5)
          annotation (Placement(transformation(extent={{-108,-12},{-88,8}})));
        Components.Controls.sens_PT sens_PT(redeclare package medium =
              ThermoCycle.Media.Toluene_CP)
          annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
        Components.Controls.sens_PT sens_PT1(redeclare package medium =
              ThermoCycle.Media.Toluene_CP)
          annotation (Placement(transformation(extent={{10,-16},{30,4}})));
        Components.source_sink_ports.OutFlow outFlow3(
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          m0=5,
          p0(displayUnit="kPa") = 696140)
          annotation (Placement(transformation(extent={{32,-4},{52,16}})));
      equation
        connect(outFlow_liquid.A,hx1counter. OutFlow_liquid) annotation (Line(
              points={{-62,46},{-74,46},{-74,16},{-48,16}}, color={28,108,200}));
        connect(inFlow_liquid.A, hx1counter.InFlow_liquid) annotation (Line(
              points={{53.4,39.8},{66.7,39.8},{66.7,16.5},{-0.48,16.5}}, color=
                {28,108,200}));
        connect(inFlow.node, sens_PT.inFlow) annotation (Line(points={{-105.6,
                -2.4},{-92.8,-2.4},{-92.8,0.2},{-80,0.2}}, color={0,0,255}));
        connect(sens_PT.outFlow, hx1counter.InFlow_working) annotation (Line(
              points={{-60,0.2},{-54,0.2},{-54,-3.5},{-48,-3.5}}, color={0,0,
                255}));
        connect(hx1counter.OutFlow_working, sens_PT1.inFlow) annotation (Line(
              points={{-0.48,-3},{4.76,-3},{4.76,-5.8},{10,-5.8}}, color={0,0,
                255}));
        connect(outFlow3.node, sens_PT1.outFlow) annotation (Line(points={{50.2,
                5.8},{57.1,5.8},{57.1,-5.8},{30,-5.8}}, color={0,0,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end testModel2;

      model testModel3
        Components.main_components.Hx1counter hx1counter(
          N=20,
          V_working=0.1753 + 1.1414,
          V_liquid=0.1224 + 0.2279,
          M=0.4*1000,
          Cp=4.2e3,
          Cp_wall=390,
          M_wall=50,
          Ul=2350,
          Uv=1700,
          redeclare package Medium = ThermoCycle.Media.Toluene_CP,
          F_working=100.05,
          F_liquid=120,
          m_working=5,
          m_liquid=25,
          pstart(displayUnit="kPa") = 696000,
          T1_liquid(displayUnit="K") = 527.91,
          Tn_liquid(displayUnit="K") = 463.65,
          Utp=8000,
          U=3000,
          Tin_start(displayUnit="K") = 383,
          Tout_start(displayUnit="K") = 470)
          annotation (Placement(transformation(extent={{-32,-30},{16,20}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid(T0(displayUnit=
               "K") = 527)
          annotation (Placement(transformation(extent={{50,20},{70,40}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid(m0=15.91)
          annotation (Placement(transformation(extent={{-46,26},{-26,46}})));
        Components.source_sink_ports.InFlow inFlow(
          UseT=true,
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          p0(displayUnit="kPa") = 696000,
          T(displayUnit="K") = 389,
          m0=-5)
          annotation (Placement(transformation(extent={{-92,-22},{-72,-2}})));
        Components.Controls.sens_PT sens_PT(redeclare package medium =
              ThermoCycle.Media.Toluene_CP) annotation (Placement(
              transformation(extent={{-64,-20},{-44,0}})));
        Components.Controls.sens_PT sens_PT1(redeclare package medium =
              ThermoCycle.Media.Toluene_CP)
          annotation (Placement(transformation(extent={{26,-26},{46,-6}})));
        Components.main_components.Turbine_simple turbine_simple1(
          n_turbine=0.681,
          n_ele=0.89,
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          pin_start(displayUnit="kPa") = 696000,
          pout_start(displayUnit="kPa") = 102000,
          Tin_start(displayUnit="K") = 468.99,
          Np=50,
          V_s=0.0063)
          annotation (Placement(transformation(extent={{48,-76},{92,-32}})));
        Components.source_sink_ports.OutFlow outFlow2(
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          m0=5,
          p0(displayUnit="kPa") = 102500)
          annotation (Placement(transformation(extent={{44,-90},{64,-70}})));
        Components.source_sink_ports.InFlow inFlow1(
          UseT=true,
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          p0(displayUnit="kPa") = 696000,
          m0=-5,
          T(displayUnit="degC") = 523.15)
          annotation (Placement(transformation(extent={{66,-24},{86,-4}})));
        Components.source_sink_ports.OutFlow outFlow1(
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          m0=5,
          p0(displayUnit="kPa") = 700000)
          annotation (Placement(transformation(extent={{16,-42},{36,-22}})));
        Components.main_components.pump_type3 pump_type3_1(
          n_pump=0.681,
          p0(displayUnit="kPa") = 249226,
          pstart_out(displayUnit="kPa") = 1440000,
          v_s=1.7423e-4,
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          T0=310.15)
          annotation (Placement(transformation(extent={{-42,-84},{-82,-34}})));
        Modelica.Blocks.Sources.Constant const(k=50)
          annotation (Placement(transformation(extent={{-102,-68},{-82,-48}})));
        Components.source_sink_ports.InFlow inFlow2(
          UseT=true,
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          m0=-5,
          p0(displayUnit="kPa") = 100000,
          T(displayUnit="degC") = 389)
          annotation (Placement(transformation(extent={{-50,-86},{-30,-66}})));
        Components.source_sink_ports.OutFlow outFlow3(
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          m0=5,
          p0(displayUnit="kPa") = 700000)
          annotation (Placement(transformation(extent={{-98,-42},{-78,-22}})));
      equation
        connect(outFlow_liquid.A,hx1counter. OutFlow_liquid) annotation (Line(
              points={{-46,36},{-70,36},{-70,6},{-32,6}},   color={28,108,200}));
        connect(inFlow_liquid.A,hx1counter. InFlow_liquid) annotation (Line(
              points={{69.4,29.8},{82.7,29.8},{82.7,6.5},{15.52,6.5}},
              color={28,108,200}));
        connect(inFlow.node,sens_PT. inFlow) annotation (Line(points={{-89.6,
                -12.4},{-76.8,-12.4},{-76.8,-9.8},{-64,-9.8}},      color={0,0,
                255}));
        connect(sens_PT.outFlow,hx1counter. InFlow_working) annotation (Line(
              points={{-44,-9.8},{-38,-9.8},{-38,-13.5},{-32,-13.5}},
              color={0,0,255}));
        connect(hx1counter.OutFlow_working,sens_PT1. inFlow) annotation (Line(
              points={{15.52,-13},{20.76,-13},{20.76,-15.8},{26,-15.8}},
              color={0,0,255}));
        connect(turbine_simple1.OutFlow,outFlow2. node) annotation (Line(points={{80.12,
                -68.52},{72.06,-68.52},{72.06,-80.2},{62.2,-80.2}},    color={0,
                0,255}));
        connect(inFlow1.node, turbine_simple1.InFlow) annotation (Line(points={
                {68.4,-14.4},{68.4,-30.2},{59,-30.2},{59,-44.76}}, color={0,0,
                255}));
        connect(outFlow1.node, sens_PT1.outFlow) annotation (Line(points={{34.2,
                -32.2},{34.2,-24.1},{46,-24.1},{46,-15.8}}, color={0,0,255}));
        connect(const.y, pump_type3_1.u) annotation (Line(points={{-81,-58},{
                -75.5,-58},{-75.5,-56},{-68,-56}}, color={0,0,127}));
        connect(pump_type3_1.InFlow, inFlow2.node) annotation (Line(points={{
                -53.6,-65},{-53.6,-70.5},{-47.6,-70.5},{-47.6,-76.4}}, color={0,
                0,255}));
        connect(outFlow3.node, pump_type3_1.OutFlow) annotation (Line(points={{
                -79.8,-32.2},{-69.8,-32.2},{-69.8,-47.5},{-68.4,-47.5}}, color=
                {0,0,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end testModel3;

      model testModel_001
        Components.main_components.Pump_liquid pump_liquid
          annotation (Placement(transformation(extent={{64,4},{44,24}})));
        Modelica.Blocks.Sources.Constant const(k=35)
          annotation (Placement(transformation(extent={{64,24},{54,34}})));
        Modelica.Blocks.Sources.Constant const1(k=800)
          annotation (Placement(transformation(extent={{-48,76},{-38,86}})));
        Modelica.Blocks.Sources.Constant const2(k=25)
          annotation (Placement(transformation(extent={{16,82},{26,92}})));
        Components.main_components.SolarCollector solarCollector1(
          K=1,
          N=20,
          M=137,
          m_flow=4.5,
          Cp=2.34e3,
          Area=3000,
          T0=403.15)
          annotation (Placement(transformation(extent={{-38,34},{0,74}})));
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Ul=1000,
          Uv=700,
          Cp_wall=0.39e3,
          max_drhode=200,
          witdth_x=0.1,
          Utp=5000,
          F_working=8.8,
          F_liquid=8.8,
          N=20,
          M_wall=10,
          M=300,
          Cp=2.34e3,
          Tin_start=312.15,
          Tout_start=388.15,
          T1_liquid=403.15,
          Tn_liquid=373.15,
          m_working=3,
          m_liquid=45)
          annotation (Placement(transformation(extent={{-26,-2},{16,32}})));
        Components.main_components.Tank tank(redeclare package medium =
              ThermoCycle.Media.R245fa_CP, T0(displayUnit="degC") = 313.15)
          annotation (Placement(transformation(extent={{-76,-72},{-34,-36}})));
        Components.main_components.Pump_type2 pump_type2_1(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          Np=35,
          T0=311.15,
          v_s=1.2801e-4)
          annotation (Placement(transformation(extent={{-106,-36},{-40,22}})));
        Components.main_components.Hx1counter condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          N=10,
          Ul=1000,
          Uv=700,
          M=700,
          Cp=4.2e3,
          M_wall=50,
          Cp_wall=0.39e3,
          witdth_x=0.1,
          Utp=5000,
          F_working=13.8,
          F_liquid=13.8,
          max_drhode=200,
          Tin_start=371.15,
          Tout_start=303.15,
          T1_liquid=303.15,
          Tn_liquid=313.15,
          m_working=3,
          m_liquid=50,
          pstart=300000)
          annotation (Placement(transformation(extent={{40,-46},{-8,-82}})));
        Components.main_components.Turbine_simple turbine_simple(redeclare
            package medium = ThermoCycle.Media.R245fa_CP,
          n_turbine=0.7,
          Np=38,
          V_s=3.121e-3,
          Tin_start=383.15)
          annotation (Placement(transformation(extent={{44,-52},{64,-32}})));
        Components.source_sink_ports.OutFlow_Tmflow outFlow_Tmflow(m0=1, use_m0=
             false,
          N=10)
          annotation (Placement(transformation(extent={{74,-96},{94,-76}})));
        Components.source_sink_ports.InFlow_Tmflow inFlow_Tmflow(
          use_m0=true,
          T0=298.15,
          m0=-120)
          annotation (Placement(transformation(extent={{-88,-92},{-68,-72}})));
        Components.main_components.sens_liquid sens_liquid
          annotation (Placement(transformation(extent={{38,44},{58,64}})));
        Components.main_components.Tank_liquid tank_liquid2(
          rou=825,
          Cp=2.34e3,
          Area=0.1,
          L0=0.1,
          T0=403.15)
          annotation (Placement(transformation(extent={{72,-14},{92,6}})));
        Components.Controls.SuperHeat superHeat(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{22,-24},{42,-4}})));
      equation
        connect(const.y,pump_liquid. u) annotation (Line(points={{53.5,29},{
                38.5,29},{38.5,17.9},{53.9,17.9}},
                                               color={0,0,127}));
        connect(const2.y,solarCollector1. Tam_in) annotation (Line(points={{26.5,87},
                {26.5,77.5},{-10.26,77.5},{-10.26,65.6}},         color={0,0,
                127}));
        connect(pump_type2_1.InFlow,tank. OutFlow) annotation (Line(points={{-86.86,
                -13.96},{-91.43,-13.96},{-91.43,-63},{-66.34,-63}},      color=
                {0,0,255}));
        connect(tank.InFlow,condenser. OutFlow_working) annotation (Line(points={{-42.4,
                -46.8},{-40,-46.8},{-40,-58.24},{-7.52,-58.24}},   color={0,0,
                255}));
        connect(turbine_simple.OutFlow,condenser. InFlow_working) annotation (
            Line(points={{58.6,-48.6},{58.6,-57.3},{40,-57.3},{40,-57.88}},
              color={0,0,255}));
        connect(hx1counter.InFlow_working,pump_type2_1. OutFlow) annotation (
            Line(points={{-26,9.22},{-48,9.22},{-48,6.34},{-62.44,6.34}}, color=
               {0,0,255}));
        connect(condenser.OutFlow_liquid,outFlow_Tmflow. inFlow) annotation (
            Line(points={{40,-71.92},{63,-71.92},{63,-86.2},{74,-86.2}}, color=
                {28,108,200}));
        connect(inFlow_Tmflow.outFlow,condenser. InFlow_liquid) annotation (
            Line(points={{-68,-82},{-38,-82},{-38,-72.28},{-7.52,-72.28}},
              color={28,108,200}));
        connect(solarCollector1.outFlow,sens_liquid. inFlow) annotation (Line(
              points={{0,54.4},{34,54.4},{34,54.2},{42.4,54.2}},  color={28,108,
                200}));
        connect(hx1counter.InFlow_liquid,pump_liquid. outFlow) annotation (Line(
              points={{15.58,22.82},{29.79,22.82},{29.79,14},{44,14}}, color={
                28,108,200}));
        connect(pump_liquid.inFlow,tank_liquid2. outFlow) annotation (Line(
              points={{64,14},{70,14},{70,-9.6},{77.2,-9.6}},
                                                            color={28,108,200}));
        connect(sens_liquid.outFlow,tank_liquid2. inFlow) annotation (Line(
              points={{53.2,54.2},{53.2,52.1},{88.2,52.1},{88.2,0.2}},  color={
                28,108,200}));
        connect(hx1counter.OutFlow_liquid,solarCollector1. inFlow) annotation (
            Line(points={{-26,22.48},{-82,22.48},{-82,52.8801},{-38,52.8801},{
                -38,54}},                                               color={
                28,108,200}));
        connect(hx1counter.OutFlow_working,superHeat. inFlow) annotation (Line(
              points={{15.58,9.56},{15.58,-3.22},{22,-3.22},{22,-13.8}},
                                                                       color={0,
                0,255}));
        connect(superHeat.outFlow,turbine_simple. InFlow) annotation (Line(
              points={{42,-13.8},{46,-13.8},{46,-37.8},{49,-37.8}},
                                                                  color={0,0,
                255}));
        connect(const1.y,solarCollector1. I_in) annotation (Line(points={{-37.5,
                81},{-37.5,79.5},{-21.28,79.5},{-21.28,65.6}},
                                                           color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end testModel_001;

      model Condenser_toluene
        Components.main_components.Hx1counter hx1counter1(
          N=20,
          V_working=0.1753 + 1.1414,
          V_liquid=0.1224 + 0.2279,
          M=0.4*1000,
          Cp=4.2e3,
          Cp_wall=390,
          M_wall=50,
          Ul=2350,
          Uv=1700,
          redeclare package Medium = ThermoCycle.Media.Toluene_CP,
          F_working=100.05,
          F_liquid=120,
          m_working=5,
          m_liquid=25,
          pstart(displayUnit="kPa") = 696000,
          T1_liquid(displayUnit="K") = 527.91,
          Tn_liquid(displayUnit="K") = 463.65,
          Utp=8000,
          U=3000,
          Tin_start(displayUnit="K") = 383,
          Tout_start(displayUnit="K") = 470)
          annotation (Placement(transformation(extent={{-42,0},{6,50}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(
                                                                   m0=15.91)
          annotation (Placement(transformation(extent={{-144,64},{-124,84}})));
        Components.Controls.sens_PT sens_PT(redeclare package medium =
              ThermoCycle.Media.Toluene_CP)
          annotation (Placement(transformation(extent={{-94,10},{-74,30}})));
        Components.Controls.sens_PT sens_PT1(redeclare package medium =
              ThermoCycle.Media.Toluene_CP)
          annotation (Placement(transformation(extent={{32,-2},{52,18}})));
        Components.main_components.Turbine_simple turbine_simple1(
          n_turbine=0.681,
          n_ele=0.89,
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          pin_start(displayUnit="kPa") = 696000,
          pout_start(displayUnit="kPa") = 102000,
          Tin_start(displayUnit="K") = 468.99,
          V_s=0.0068,
          Np=45)
          annotation (Placement(transformation(extent={{66,-16},{110,28}})));
        Components.main_components.SolarCollector solarCollector(
          N=20,
          K=1,
          Cp=2400,
          m_flow=25,
          M=300,
          Area=5900,
          T0=473.15)
          annotation (Placement(transformation(extent={{-30,64},{32,110}})));
        Modelica.Blocks.Sources.Constant const1(k=600)
          annotation (Placement(transformation(extent={{-98,106},{-78,126}})));
        Modelica.Blocks.Sources.Constant const2(k=25)
          annotation (Placement(transformation(extent={{12,124},{32,144}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid2(
                                                                   m0=15.91)
          annotation (Placement(transformation(extent={{96,62},{116,82}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid2(
                                                                 T0(displayUnit=
               "K") = 527)
          annotation (Placement(transformation(extent={{58,34},{78,54}})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT
          annotation (Placement(transformation(extent={{-70,80},{-50,100}})));
        Components.main_components.sens_liquid sens_liquid
          annotation (Placement(transformation(extent={{-82,66},{-102,86}})));
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.Toluene_CP,
          M=700,
          Cp=2400,
          M_wall=50,
          Cp_wall=380,
          m_working=5,
          m_liquid=20,
          V_working=1,
          V_liquid=1,
          U=1100,
          N=20,
          max_drhode=10,
          F_working=160,
          F_liquid=160,
          witdth_x=0.05,
          Mcons=true,
          Tin_start=453.15,
          Tout_start=373.15,
          T1_liquid=343.15,
          Tn_liquid=403.15,
          pstart=110000)
          annotation (Placement(transformation(extent={{40,-6},{-36,-60}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid1(T0(
              displayUnit="degC") = 343.15)
          annotation (Placement(transformation(extent={{-66,-76},{-46,-56}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid3(m0=30)
          annotation (Placement(transformation(extent={{66,-54},{86,-34}})));
        Components.main_components.pump_type3 pump_type3_1(
          n_pump=0.681,
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          p0(displayUnit="kPa") = 102500,
          pstart_out(displayUnit="kPa") = 696140,
          T0(displayUnit="K") = 383,
          v_s=1.2867e-4)
          annotation (Placement(transformation(extent={{-74,-36},{-114,14}})));
        Modelica.Blocks.Sources.Constant const(k=45)
          annotation (Placement(transformation(extent={{-164,-20},{-144,0}})));
        Components.main_components.Tank tank(
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          Aera=1,
          L0=1,
          p0=110000,
          T0=353.15)
          annotation (Placement(transformation(extent={{-72,-42},{-38,-12}})));
        Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(
          tableOnFile=true,
          columns={2},
          tableName="nana",
          fileName="nana.mat")
          annotation (Placement(transformation(extent={{-120,144},{-106,158}})));
        Modelica.Blocks.Math.Gain gain(k=0.8)
          annotation (Placement(transformation(extent={{-90,138},{-70,158}})));
      equation
        connect(sens_PT.outFlow, hx1counter1.InFlow_working) annotation (Line(
              points={{-74,20.2},{-60,20.2},{-60,16.5},{-42,16.5}}, color={0,0,
                255}));
        connect(hx1counter1.OutFlow_working, sens_PT1.inFlow) annotation (Line(
              points={{5.52,17},{10.76,17},{10.76,8.2},{32,8.2}}, color={0,0,
                255}));
        connect(sens_PT1.outFlow, turbine_simple1.InFlow) annotation (Line(
              points={{52,8.2},{70,8.2},{70,15.24},{77,15.24}}, color={0,0,255}));
        connect(const2.y, solarCollector.Tam_in) annotation (Line(points={{33,
                134},{42,134},{42,100.34},{15.26,100.34}}, color={0,0,127}));
        connect(outFlow_liquid2.A, inFlow_liquid2.A) annotation (Line(points={{
                96,72},{88,72},{88,43.8},{77.4,43.8}}, color={28,108,200}));
        connect(solarCollector.outFlow, hx1counter1.InFlow_liquid) annotation (
            Line(points={{32,87.46},{44,87.46},{44,36.5},{5.52,36.5}}, color={
                28,108,200}));
        connect(inFlow_liquid_useT.A, solarCollector.inFlow) annotation (Line(
              points={{-50.6,89.8},{-45.3,89.8},{-45.3,87},{-30,87}}, color={28,
                108,200}));
        connect(sens_liquid.inFlow, hx1counter1.OutFlow_liquid) annotation (
            Line(points={{-86.4,76.2},{-86.4,58.1},{-42,58.1},{-42,36}}, color=
                {28,108,200}));
        connect(sens_liquid.outFlow, outFlow_liquid1.A) annotation (Line(points=
               {{-97.2,76.2},{-111.4,76.2},{-111.4,74},{-144,74}}, color={28,
                108,200}));
        connect(sens_liquid.y, inFlow_liquid_useT.u) annotation (Line(points={{-91.6,
                79.2},{-91.8,79.2},{-91.8,94},{-61.2,94}},         color={0,0,
                127}));
        connect(hx1counter.InFlow_working, turbine_simple1.OutFlow) annotation (
           Line(points={{40,-23.82},{68,-23.82},{68,-8.52},{98.12,-8.52}},
              color={0,0,255}));
        connect(inFlow_liquid1.A, hx1counter.InFlow_liquid) annotation (Line(
              points={{-46.6,-66.2},{-49.3,-66.2},{-49.3,-45.42},{-35.24,-45.42}},
              color={28,108,200}));
        connect(hx1counter.OutFlow_liquid, outFlow_liquid3.A) annotation (Line(
              points={{40,-44.88},{54,-44.88},{54,-44},{66,-44}}, color={28,108,
                200}));
        connect(const.y, pump_type3_1.u) annotation (Line(points={{-143,-10},{
                -100,-10},{-100,-8}}, color={0,0,127}));
        connect(tank.InFlow, hx1counter.OutFlow_working) annotation (Line(
              points={{-44.8,-21},{-44.8,-18.5},{-35.24,-18.5},{-35.24,-24.36}},
              color={0,0,255}));
        connect(tank.OutFlow, pump_type3_1.InFlow) annotation (Line(points={{-64.18,
                -34.5},{-72.09,-34.5},{-72.09,-17},{-85.6,-17}},        color={
                0,0,255}));
        connect(pump_type3_1.OutFlow, sens_PT.inFlow) annotation (Line(points={
                {-100.4,0.5},{-100.4,13.25},{-94,13.25},{-94,20.2}}, color={0,0,
                255}));
        connect(combiTimeTable.y[1], gain.u) annotation (Line(points={{-105.3,
                151},{-98.65,151},{-98.65,148},{-92,148}}, color={0,0,127}));
        connect(gain.y, solarCollector.I_in) annotation (Line(points={{-69,148},
                {-4,148},{-4,100.34},{-2.72,100.34}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{58,106},{96,70}},
                lineColor={28,108,200},
                textString="注意：要用最后一个求解算法
并且把离散间隔取尽量小
误差取大一些才能计算")}));
      end Condenser_toluene;

      model TT
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          N=20,
          M=500,
          Cp=2400,
          M_wall=50,
          Cp_wall=380,
          Tin_start=343.15,
          Tout_start=383.15,
          T1_liquid=423.15,
          Tn_liquid=353.15,
          m_working=2,
          m_liquid=4,
          pstart(displayUnit="kPa") = 700000)
          annotation (Placement(transformation(extent={{-42,-8},{54,46}})));
        Components.source_sink_ports.InFlow inFlow(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0(displayUnit="kPa") = 700000,
          UseT=true,
          m0=-3,
          T=343.15)
          annotation (Placement(transformation(extent={{-48,-26},{-28,-6}})));
        Components.source_sink_ports.OutFlow outFlow(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          m0=3,
          p0(displayUnit="kPa") = 200000)
          annotation (Placement(transformation(extent={{24,-92},{44,-72}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid(T0=423.15)
          annotation (Placement(transformation(extent={{34,54},{54,74}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid(m0=5)
          annotation (Placement(transformation(extent={{-40,54},{-20,74}})));
        Components.main_components.Turbine_simple turbine_simple(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          Np=55,
          V_s=0.0007)
          annotation (Placement(transformation(extent={{62,-76},{116,6}})));
        Modelica.Blocks.Sources.CombiTimeTable combiTimeTable1(
          tableOnFile=true,
          columns={2},
          tableName="nana",
          fileName="nana.mat")
          annotation (Placement(transformation(extent={{-92,48},{-78,62}})));
      equation
        connect(inFlow.node, hx1counter.InFlow_working) annotation (Line(points=
               {{-45.6,-16.4},{-45.6,-3.2},{-42,-3.2},{-42,9.82}}, color={0,0,
                255}));
        connect(inFlow_liquid.A, hx1counter.InFlow_liquid) annotation (Line(
              points={{53.4,63.8},{53.4,47.9},{53.04,47.9},{53.04,31.42}},
              color={28,108,200}));
        connect(outFlow_liquid.A, hx1counter.OutFlow_liquid) annotation (Line(
              points={{-40,64},{-44,64},{-44,30.88},{-42,30.88}}, color={28,108,
                200}));
        connect(hx1counter.OutFlow_working, turbine_simple.InFlow) annotation (
            Line(points={{53.04,10.36},{53.04,-3.82},{75.5,-3.82},{75.5,-17.78}},
              color={0,0,255}));
        connect(outFlow.node, turbine_simple.OutFlow) annotation (Line(points={
                {42.2,-82.2},{72.1,-82.2},{72.1,-62.06},{101.42,-62.06}}, color=
               {0,0,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TT;

      model solar_tolune_final
        Components.main_components.Hx1counter hx1counter1(
          N=20,
          V_working=0.1753 + 1.1414,
          V_liquid=0.1224 + 0.2279,
          M=0.4*1000,
          Cp=4.2e3,
          Cp_wall=390,
          M_wall=50,
          Uv=1700,
          redeclare package Medium = ThermoCycle.Media.Toluene_CP,
          m_working=5,
          m_liquid=25,
          pstart(displayUnit="kPa") = 696000,
          T1_liquid(displayUnit="K") = 527.91,
          Tn_liquid(displayUnit="K") = 463.65,
          Tin_start(displayUnit="K") = 383,
          Tout_start(displayUnit="K") = 470,
          F_working=130.61,
          F_liquid=130.61,
          U=2000,
          Ul=1050,
          Utp=4000)
          annotation (Placement(transformation(extent={{-24,-44},{20,2}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(
                                                                   m0=15.91)
          annotation (Placement(transformation(extent={{-112,-6},{-92,14}})));
        Components.Controls.sens_PT sens_PT(redeclare package medium =
              ThermoCycle.Media.Toluene_CP)
          annotation (Placement(transformation(extent={{-88,-30},{-68,-10}})));
        Components.Controls.sens_PT sens_PT1(redeclare package medium =
              ThermoCycle.Media.Toluene_CP)
          annotation (Placement(transformation(extent={{40,-40},{60,-20}})));
        Components.main_components.Turbine_simple turbine_simple1(
          n_turbine=0.681,
          n_ele=0.89,
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          pin_start(displayUnit="kPa") = 696000,
          Tin_start(displayUnit="K") = 468.99,
          V_s=0.0068,
          Np=45,
          pout_start(displayUnit="kPa") = 110000)
          annotation (Placement(transformation(extent={{72,-56},{116,-12}})));
        Modelica.Blocks.Sources.Constant const2(k=25)
          annotation (Placement(transformation(extent={{92,80},{112,100}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid2(
                                                                   m0=15.91)
          annotation (Placement(transformation(extent={{102,22},{122,42}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid2(
                                                                 T0(displayUnit=
               "K") = 527)
          annotation (Placement(transformation(extent={{66,22},{86,42}})));
        Components.main_components.sens_liquid sens_liquid
          annotation (Placement(transformation(extent={{-60,-16},{-80,4}})));
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.Toluene_CP,
          M=700,
          Cp=2400,
          M_wall=50,
          Cp_wall=380,
          m_working=5,
          m_liquid=20,
          V_working=1,
          V_liquid=1,
          N=20,
          F_working=174.73,
          F_liquid=174.73,
          Mcons=true,
          max_drhode=100,
          witdth_x=0.05,
          Ul=1000,
          Uv=2000,
          U=1100,
          Tin_start=453.15,
          Tout_start=373.15,
          T1_liquid=343.15,
          Tn_liquid=403.15,
          pstart=110000)
          annotation (Placement(transformation(extent={{28,-62},{-30,-100}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid1(T0(
              displayUnit="degC") = 343.15)
          annotation (Placement(transformation(extent={{-78,-98},{-58,-78}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid3(m0=30)
          annotation (Placement(transformation(extent={{72,-94},{92,-74}})));
        Components.main_components.pump_type3 pump_type3_1(
          n_pump=0.681,
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          p0(displayUnit="kPa") = 102500,
          pstart_out(displayUnit="kPa") = 696140,
          T0(displayUnit="K") = 383,
          v_s=1.2867e-4)
          annotation (Placement(transformation(extent={{-68,-76},{-108,-26}})));
        Modelica.Blocks.Sources.Constant const(k=51)
          annotation (Placement(transformation(extent={{-158,-60},{-138,-40}})));
        Components.main_components.Tank tank(
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          Aera=1,
          L0=1,
          p0=110000,
          T0=353.15)
          annotation (Placement(transformation(extent={{-72,-78},{-38,-48}})));
        Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(
          tableOnFile=true,
          columns={2},
          tableName="nana",
          fileName="nana.mat")
          annotation (Placement(transformation(extent={{-114,98},{-100,112}})));
        Modelica.Blocks.Math.Gain gain(k=0.8)
          annotation (Placement(transformation(extent={{-64,92},{-44,112}})));
        Components.main_components.SolarCollector_New solarCollector_New(
          N=20,
          K=1,
          M=300,
          Cp=2400,
          m_flow=25,
          Date=180,
          Area=6000,
          T0=473.15)
          annotation (Placement(transformation(extent={{-22,16},{24,54}})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT1
          annotation (Placement(transformation(extent={{-62,20},{-42,40}})));
        Modelica.Blocks.Sources.Constant const1(k=450)
          annotation (Placement(transformation(extent={{-84,56},{-64,76}})));
      equation
        connect(sens_PT.outFlow,hx1counter1. InFlow_working) annotation (Line(
              points={{-68,-19.8},{-54,-19.8},{-54,-28.82},{-24,-28.82}},
                                                                    color={0,0,
                255}));
        connect(hx1counter1.OutFlow_working,sens_PT1. inFlow) annotation (Line(
              points={{19.56,-28.36},{28.76,-28.36},{28.76,-29.8},{40,-29.8}},
                                                                  color={0,0,
                255}));
        connect(sens_PT1.outFlow,turbine_simple1. InFlow) annotation (Line(
              points={{60,-29.8},{76,-29.8},{76,-24.76},{83,-24.76}},
                                                                color={0,0,255}));
        connect(outFlow_liquid2.A,inFlow_liquid2. A) annotation (Line(points={{102,32},
                {94,32},{94,31.8},{85.4,31.8}},        color={28,108,200}));
        connect(sens_liquid.inFlow,hx1counter1. OutFlow_liquid) annotation (
            Line(points={{-64.4,-5.8},{-64.4,-7.9},{-24,-7.9},{-24,-10.88}},
                                                                         color=
                {28,108,200}));
        connect(sens_liquid.outFlow,outFlow_liquid1. A) annotation (Line(points={{-75.2,
                -5.8},{-83.4,-5.8},{-83.4,4},{-112,4}},            color={28,
                108,200}));
        connect(hx1counter.InFlow_working,turbine_simple1. OutFlow) annotation (
           Line(points={{28,-74.54},{74,-74.54},{74,-48.52},{104.12,-48.52}},
              color={0,0,255}));
        connect(inFlow_liquid1.A,hx1counter. InFlow_liquid) annotation (Line(
              points={{-58.6,-88.2},{-43.3,-88.2},{-43.3,-89.74},{-29.42,-89.74}},
              color={28,108,200}));
        connect(hx1counter.OutFlow_liquid,outFlow_liquid3. A) annotation (Line(
              points={{28,-89.36},{60,-89.36},{60,-84},{72,-84}}, color={28,108,
                200}));
        connect(const.y,pump_type3_1. u) annotation (Line(points={{-137,-50},{
                -94,-50},{-94,-48}},  color={0,0,127}));
        connect(tank.InFlow,hx1counter. OutFlow_working) annotation (Line(
              points={{-44.8,-57},{-44.8,-58.5},{-29.42,-58.5},{-29.42,-74.92}},
              color={0,0,255}));
        connect(tank.OutFlow,pump_type3_1. InFlow) annotation (Line(points={{-64.18,
                -70.5},{-66.09,-70.5},{-66.09,-57},{-79.6,-57}},        color={
                0,0,255}));
        connect(pump_type3_1.OutFlow,sens_PT. inFlow) annotation (Line(points={{-94.4,
                -39.5},{-94.4,-26.75},{-88,-26.75},{-88,-19.8}},     color={0,0,
                255}));
        connect(combiTimeTable.y[1],gain. u) annotation (Line(points={{-99.3,
                105},{-84.65,105},{-84.65,102},{-66,102}}, color={0,0,127}));
        connect(const2.y, solarCollector_New.Tam_in) annotation (Line(points={{113,90},
                {10,90},{10,46.02},{11.58,46.02}},          color={0,0,127}));
        connect(inFlow_liquid_useT1.A, solarCollector_New.inFlow) annotation (
            Line(points={{-42.6,29.8},{-36.3,29.8},{-36.3,35},{-22,35}}, color=
                {28,108,200}));
        connect(sens_liquid.y, inFlow_liquid_useT1.u) annotation (Line(points={{-69.6,
                -2.8},{-69.6,35.6},{-53.2,35.6},{-53.2,34}},           color={0,
                0,127}));
        connect(gain.y, solarCollector_New.I_in) annotation (Line(points={{-43,102},
                {-2,102},{-2,46.02},{-1.76,46.02}},   color={0,0,127}));
        connect(solarCollector_New.outFlow, hx1counter1.InFlow_liquid)
          annotation (Line(points={{24,35.38},{48,35.38},{48,-10.42},{19.56,
                -10.42}}, color={28,108,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{64,66},{102,30}},
                lineColor={28,108,200},
                textString="注意：要用最后一个求解算法
并且把离散间隔取尽量小
误差取大一些才能计算")}));
      end solar_tolune_final;

      model SORC_toluene
        Components.main_components.Hx1counter hx1counter2(
          N=20,
          V_working=0.1753 + 1.1414,
          V_liquid=0.1224 + 0.2279,
          M=0.4*1000,
          Cp=4.2e3,
          Cp_wall=390,
          M_wall=50,
          Ul=2350,
          Uv=1700,
          redeclare package Medium = ThermoCycle.Media.Toluene_CP,
          F_working=100.05,
          F_liquid=120,
          m_working=5,
          m_liquid=25,
          T1_liquid(displayUnit="K") = 527.91,
          Tn_liquid(displayUnit="K") = 463.65,
          Utp=8000,
          U=3000,
          Tin_start(displayUnit="K") = 383,
          Tout_start(displayUnit="K") = 470,
          pstart(displayUnit="kPa") = 696000)
          annotation (Placement(transformation(extent={{-16,-26},{32,24}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid4(
                                                                   m0=15.91)
          annotation (Placement(transformation(extent={{-72,6},{-52,26}})));
        Components.Controls.sens_PT sens_PT2(
                                            redeclare package medium =
              ThermoCycle.Media.Toluene_CP)
          annotation (Placement(transformation(extent={{-68,-16},{-48,4}})));
        Components.Controls.sens_PT sens_PT3(redeclare package medium =
              ThermoCycle.Media.Toluene_CP)
          annotation (Placement(transformation(extent={{58,-28},{78,-8}})));
        Components.main_components.Turbine_simple turbine_simple2(
          n_turbine=0.681,
          n_ele=0.89,
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          pin_start(displayUnit="kPa") = 696000,
          pout_start(displayUnit="kPa") = 102000,
          Tin_start(displayUnit="K") = 468.99,
          V_s=0.0058,
          Np=45)
          annotation (Placement(transformation(extent={{92,-42},{136,2}})));
        Components.main_components.SolarCollector solarCollector(
          N=20,
          K=1,
          Cp=2400,
          m_flow=25,
          M=300,
          Area=5900,
          T0=473.15)
          annotation (Placement(transformation(extent={{-4,38},{58,84}})));
        Modelica.Blocks.Sources.Constant const3(k=600)
          annotation (Placement(transformation(extent={{-94,66},{-74,86}})));
        Modelica.Blocks.Sources.Constant const4(k=25)
          annotation (Placement(transformation(extent={{28,184},{48,204}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid5(
                                                                   m0=15.91)
          annotation (Placement(transformation(extent={{122,36},{142,56}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid3(
                                                                 T0(displayUnit=
               "K") = 527)
          annotation (Placement(transformation(extent={{84,8},{104,28}})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT1
          annotation (Placement(transformation(extent={{-42,48},{-22,68}})));
        Components.main_components.sens_liquid sens_liquid1
          annotation (Placement(transformation(extent={{-56,40},{-76,60}})));
        Components.main_components.pump_type3 pump_type3_2(
          n_pump=0.681,
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          p0(displayUnit="kPa") = 102500,
          pstart_out(displayUnit="kPa") = 696140,
          T0(displayUnit="K") = 383,
          v_s=1.2867e-4)
          annotation (Placement(transformation(extent={{-48,-62},{-88,-12}})));
        Modelica.Blocks.Sources.Constant const5(k=45)
          annotation (Placement(transformation(extent={{-110,-44},{-90,-24}})));
        Modelica.Blocks.Sources.CombiTimeTable combiTimeTable1(
          tableOnFile=true,
          columns={2},
          tableName="nana",
          fileName="nana.mat")
          annotation (Placement(transformation(extent={{-96,176},{-82,190}})));
        Modelica.Blocks.Math.Gain gain1(
                                       k=0.8)
          annotation (Placement(transformation(extent={{-58,174},{-38,194}})));
        Components.main_components.SolarCollector_New solarCollector_New1(
          N=20,
          K=1,
          Cp=2400,
          m_flow=25,
          ti=10,
          M=300,
          use_eff2=true,
          Area=7188.1,
          Date=100,
          T0=473.15,
          lati=37)
          annotation (Placement(transformation(extent={{-12,96},{46,140}})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT2
          annotation (Placement(transformation(extent={{-52,110},{-32,130}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(
                                                                   m0=15.91)
          annotation (Placement(transformation(extent={{130,104},{150,124}})));
        Components.main_components.Hx1counter hx1counter(
          redeclare package Medium = ThermoCycle.Media.Toluene_CP,
          M=700,
          Cp=2400,
          M_wall=50,
          Cp_wall=380,
          m_working=5,
          m_liquid=20,
          V_working=1,
          V_liquid=1,
          U=1100,
          N=20,
          max_drhode=10,
          F_working=160,
          F_liquid=160,
          witdth_x=0.05,
          Mcons=true,
          Tin_start=453.15,
          Tout_start=373.15,
          T1_liquid=343.15,
          Tn_liquid=403.15,
          pstart=110000)
          annotation (Placement(transformation(extent={{66,-54},{-10,-108}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid1(T0(
              displayUnit="degC") = 343.15)
          annotation (Placement(transformation(extent={{-40,-124},{-20,-104}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid3(m0=30)
          annotation (Placement(transformation(extent={{92,-102},{112,-82}})));
        Components.main_components.Tank tank(
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          Aera=1,
          L0=1,
          p0=110000,
          T0=353.15)
          annotation (Placement(transformation(extent={{-46,-90},{-12,-60}})));
      equation
        connect(sens_PT2.outFlow, hx1counter2.InFlow_working) annotation (Line(
              points={{-48,-5.8},{-34,-5.8},{-34,-9.5},{-16,-9.5}}, color={0,0,
                255}));
        connect(hx1counter2.OutFlow_working,sens_PT3. inFlow) annotation (Line(
              points={{31.52,-9},{36.76,-9},{36.76,-17.8},{58,-17.8}},
                                                                  color={0,0,
                255}));
        connect(sens_PT3.outFlow,turbine_simple2. InFlow) annotation (Line(
              points={{78,-17.8},{96,-17.8},{96,-10.76},{103,-10.76}},
                                                                color={0,0,255}));
        connect(const4.y,solarCollector. Tam_in) annotation (Line(points={{49,194},
                {68,194},{68,74.34},{41.26,74.34}},        color={0,0,127}));
        connect(outFlow_liquid5.A,inFlow_liquid3. A) annotation (Line(points={{122,46},
                {114,46},{114,17.8},{103.4,17.8}},     color={28,108,200}));
        connect(inFlow_liquid_useT1.A, solarCollector.inFlow) annotation (Line(
              points={{-22.6,57.8},{-19.3,57.8},{-19.3,61},{-4,61}}, color={28,
                108,200}));
        connect(sens_liquid1.inFlow, hx1counter2.OutFlow_liquid) annotation (
            Line(points={{-60.4,50.2},{-60.4,32.1},{-16,32.1},{-16,10}}, color=
                {28,108,200}));
        connect(sens_liquid1.outFlow, outFlow_liquid4.A) annotation (Line(
              points={{-71.2,50.2},{-85.4,50.2},{-85.4,16},{-72,16}}, color={28,
                108,200}));
        connect(sens_liquid1.y, inFlow_liquid_useT1.u) annotation (Line(points=
                {{-65.6,53.2},{-59.8,53.2},{-59.8,62},{-33.2,62}}, color={0,0,
                127}));
        connect(const5.y, pump_type3_2.u)
          annotation (Line(points={{-89,-34},{-74,-34}}, color={0,0,127}));
        connect(pump_type3_2.OutFlow, sens_PT2.inFlow) annotation (Line(points=
                {{-74.4,-25.5},{-74.4,-12.75},{-68,-12.75},{-68,-5.8}}, color={
                0,0,255}));
        connect(combiTimeTable1.y[1], gain1.u) annotation (Line(points={{-81.3,
                183},{-72.65,183},{-72.65,184},{-60,184}}, color={0,0,127}));
        connect(gain1.y, solarCollector.I_in) annotation (Line(points={{-37,184},
                {22,184},{22,74.34},{23.28,74.34}}, color={0,0,127}));
        connect(sens_liquid1.y, inFlow_liquid_useT2.u) annotation (Line(points=
                {{-65.6,53.2},{-65.6,89.6},{-43.2,89.6},{-43.2,124}}, color={0,
                0,127}));
        connect(inFlow_liquid_useT2.A, solarCollector_New1.inFlow) annotation (
            Line(points={{-32.6,119.8},{-26.3,119.8},{-26.3,118},{-12,118}},
              color={28,108,200}));
        connect(gain1.y, solarCollector_New1.I_in) annotation (Line(points={{
                -37,184},{4,184},{4,130.76},{13.52,130.76}}, color={0,0,127}));
        connect(const4.y, solarCollector_New1.Tam_in) annotation (Line(points={
                {49,194},{60,194},{60,130.76},{30.34,130.76}}, color={0,0,127}));
        connect(solarCollector.outFlow, outFlow_liquid1.A) annotation (Line(
              points={{58,61.46},{86,61.46},{86,114},{130,114}}, color={28,108,
                200}));
        connect(solarCollector_New1.outFlow, hx1counter2.InFlow_liquid)
          annotation (Line(points={{46,118.44},{46,98.22},{31.52,98.22},{31.52,
                10.5}}, color={28,108,200}));
        connect(inFlow_liquid1.A,hx1counter. InFlow_liquid) annotation (Line(
              points={{-20.6,-114.2},{-23.3,-114.2},{-23.3,-93.42},{-9.24,
                -93.42}},
              color={28,108,200}));
        connect(hx1counter.OutFlow_liquid,outFlow_liquid3. A) annotation (Line(
              points={{66,-92.88},{80,-92.88},{80,-92},{92,-92}}, color={28,108,
                200}));
        connect(tank.InFlow,hx1counter. OutFlow_working) annotation (Line(
              points={{-18.8,-69},{-18.8,-66.5},{-9.24,-66.5},{-9.24,-72.36}},
              color={0,0,255}));
        connect(pump_type3_2.InFlow, tank.OutFlow) annotation (Line(points={{
                -59.6,-43},{-59.6,-62.5},{-38.18,-62.5},{-38.18,-82.5}}, color=
                {0,0,255}));
        connect(turbine_simple2.OutFlow, hx1counter.InFlow_working) annotation (
           Line(points={{124.12,-34.52},{96.06,-34.52},{96.06,-71.82},{66,
                -71.82}}, color={0,0,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{66,62},{104,26}},
                lineColor={28,108,200},
                textString="不要改参数，现在可以调节通顺")}));
      end SORC_toluene;

      model SORC_fin
        Components.main_components.Hx1counter hx1counter2(
          N=20,
          V_working=0.1753 + 1.1414,
          V_liquid=0.1224 + 0.2279,
          M=0.4*1000,
          Cp_wall=390,
          M_wall=50,
          redeclare package Medium = ThermoCycle.Media.Toluene_CP,
          Cp=2764,
          Mcons=false,
          Tout_start(displayUnit="degC") = 543.15,
          m_working=4.28,
          m_liquid=13.618,
          pstart(displayUnit="kPa") = 1600000,
          witdth_x=0.1,
          Ul=1350,
          Uv=1100,
          F_working=130.61,
          F_liquid=130.61,
          T1_liquid(displayUnit="degC") = 553.15,
          Tn_liquid(displayUnit="degC") = 473.15,
          Tin_start(displayUnit="degC") = 433.15,
          Utp=3000,
          U=1000)
          annotation (Placement(transformation(extent={{-36,-18},{12,32}})));
        Components.source_sink_ports.InFlow inFlow(
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          UseT=true,
          T(displayUnit="degC") = 423.15,
          m0=-4.8,
          p0=1600000)
          annotation (Placement(transformation(extent={{-70,-34},{-50,-14}})));
        Components.source_sink_ports.OutFlow outFlow(
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          m0=4.28,
          p0(displayUnit="kPa") = 67000)
          annotation (Placement(transformation(extent={{28,-78},{48,-58}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid4(m0=13.618)
          annotation (Placement(transformation(extent={{-112,6},{-132,26}})));
        Components.main_components.SolarCollector_New solarCollector_New1(
          N=20,
          K=1,
          Cp=2400,
          m_flow=25,
          ti=10,
          M=300,
          use_eff2=true,
          Area=7188.1,
          lati=39.08,
          Date=212,
          T0=513.15)
          annotation (Placement(transformation(extent={{-26,36},{32,80}})));
        Modelica.Blocks.Sources.Constant const4(k=25)
          annotation (Placement(transformation(extent={{48,86},{56,94}})));
        Modelica.Blocks.Sources.CombiTimeTable combiTimeTable1(
          tableOnFile=true,
          columns={2},
          tableName="nana",
          fileName="nana.mat")
          annotation (Placement(transformation(extent={{-56,82},{-44,94}})));
        Modelica.Blocks.Math.Gain gain1(
                                       k=0.8)
          annotation (Placement(transformation(extent={{-22,84},{-12,94}})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT1
          annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
        Components.main_components.sens_liquid sens_liquid1
          annotation (Placement(transformation(extent={{-80,6},{-100,26}})));
        Components.main_components.Turbine_simple turbine_simple2(
          n_ele=0.89,
          redeclare package medium = ThermoCycle.Media.Toluene_CP,
          Tin_start(displayUnit="K") = 468.99,
          Np=50,
          pin_start(displayUnit="kPa") = 1600000,
          pout_start(displayUnit="kPa") = 67000,
          n_turbine=0.8,
          V_s=0.0029)
          annotation (Placement(transformation(extent={{62,-52},{106,-8}})));
        Components.Controls.sens_PT sens_PT3(redeclare package medium =
              ThermoCycle.Media.Toluene_CP)
          annotation (Placement(transformation(extent={{22,-10},{42,10}})));
        Components.Controls.SuperHeat superHeat(redeclare package medium =
              ThermoCycle.Media.Toluene_CP)
          annotation (Placement(transformation(extent={{58,-8},{78,12}})));
      equation
        connect(inFlow.node, hx1counter2.InFlow_working) annotation (Line(
              points={{-67.6,-24.4},{-50.8,-24.4},{-50.8,-1.5},{-36,-1.5}},
              color={0,0,255}));
        connect(solarCollector_New1.outFlow, hx1counter2.InFlow_liquid)
          annotation (Line(points={{32,58.44},{56,58.44},{56,18.5},{11.52,18.5}},
              color={28,108,200}));
        connect(combiTimeTable1.y[1], gain1.u) annotation (Line(points={{-43.4,
                88},{-43.4,89},{-23,89}}, color={0,0,127}));
        connect(gain1.y, solarCollector_New1.I_in) annotation (Line(points={{
                -11.5,89},{-1.75,89},{-1.75,70.76},{-0.48,70.76}}, color={0,0,
                127}));
        connect(const4.y, solarCollector_New1.Tam_in) annotation (Line(points={
                {56.4,90},{12,90},{12,70.76},{16.34,70.76}}, color={0,0,127}));
        connect(outFlow_liquid4.A, sens_liquid1.outFlow) annotation (Line(
              points={{-112,16},{-96,16},{-96,16.2},{-95.2,16.2}}, color={28,
                108,200}));
        connect(sens_liquid1.inFlow, hx1counter2.OutFlow_liquid) annotation (
            Line(points={{-84.4,16.2},{-52.2,16.2},{-52.2,18},{-36,18}}, color=
                {28,108,200}));
        connect(sens_liquid1.y, inFlow_liquid_useT1.u) annotation (Line(points=
                {{-89.6,19.2},{-89.6,58.6},{-71.2,58.6},{-71.2,54}}, color={0,0,
                127}));
        connect(inFlow_liquid_useT1.A, solarCollector_New1.inFlow) annotation (
            Line(points={{-60.6,49.8},{-46.3,49.8},{-46.3,58},{-26,58}}, color=
                {28,108,200}));
        connect(hx1counter2.OutFlow_working, sens_PT3.inFlow) annotation (Line(
              points={{11.52,-1},{21.76,-1},{21.76,0.2},{22,0.2}}, color={0,0,
                255}));
        connect(outFlow.node, turbine_simple2.OutFlow) annotation (Line(points=
                {{46.2,-68.2},{65.1,-68.2},{65.1,-44.52},{94.12,-44.52}}, color=
               {0,0,255}));
        connect(sens_PT3.outFlow, superHeat.inFlow) annotation (Line(points={{
                42,0.2},{52,0.2},{52,2.2},{58,2.2}}, color={0,0,255}));
        connect(superHeat.outFlow, turbine_simple2.InFlow) annotation (Line(
              points={{78,2.2},{76,2.2},{76,-20.76},{73,-20.76}}, color={0,0,
                255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end SORC_fin;
    end TestToluene_Cycle;

    package WHR_system
      model WHR_ORC_system
        Components.Controls.sens_PT sens_pre_out(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{-72,46},{-52,66}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(m0=2.1)
          annotation (Placement(transformation(extent={{-92,4},{-112,24}})));
        Components.main_components.Hx1counter preHeater(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Uv=800,
          Cp_wall=400,
          N=20,
          Cp=4.2e3,
          m_liquid=2.1,
          F_working=4.3,
          F_liquid=4.3,
          Utp=1300,
          M=50,
          M_wall=10,
          m_working=0.88,
          pstart(displayUnit="kPa") = 2428700,
          Ul=3200,
          U=2500,
          Tin_start=296.45,
          Tout_start=365.35,
          T1_liquid=366.15,
          Tn_liquid=358.45)
                          annotation (Placement(transformation(
              extent={{-10,-11},{10,11}},
              rotation=90,
              origin={-73,38})));
        Components.Controls.sens_PT sens_pump_out(redeclare package medium =
              ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
              extent={{6,7},{-6,-7}},
              rotation=-90,
              origin={-71,12})));
        Modelica.Blocks.Sources.Constant const2(k=50)
          annotation (Placement(transformation(extent={{-116,-8},{-106,2}})));
        Components.main_components.Tank_rou tank_rou(
          L0=0.65,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0(displayUnit="kPa") = 270020,
          Aera=1,
          T0=296.15)
                    annotation (Placement(transformation(extent={{-52,-76},{-20,
                  -46}})));
        Components.main_components.pump_type3 pump_type3_1(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0(displayUnit="kPa") = 270420,
          pstart_out(displayUnit="kPa") = 270420,
          n_pump=0.7,
          T0=310.15,
          v_s=1.3e-5)
          annotation (Placement(transformation(extent={{-66,-46},{-98,-4}})));
        Components.main_components.Hx1counter Evaporater(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Cp_wall=400,
          N=20,
          M_wall=13,
          m_working=0.88,
          m_liquid=0.298,
          M=5,
          Utp=4000,
          U=2560,
          Cp=1400,
          pstart(displayUnit="kPa") = 2000000,
          Uv=3682,
          Ul=5181,
          F_working=0.6,
          F_liquid=0.6,
          Tin_start=391.15,
          Tout_start(displayUnit="degC") = 433,
          T1_liquid=808.15,
          Tn_liquid=423.15)
          annotation (Placement(transformation(extent={{-8,46},{24,74}})));
        Components.main_components.Turbine_simple_control
          turbine_simple_control(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          pin_start(displayUnit="kPa") = 2828000,
          pout_start(displayUnit="kPa") = 270000,
          V_s=1.9e-4,
          Tin_start=453.15)
          annotation (Placement(transformation(extent={{46,-6},{76,28}})));
        Modelica.Blocks.Sources.Constant const1(k=40)
          annotation (Placement(transformation(extent={{70,40},{80,50}})));
        Components.main_components.Hx1counter Condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Cp_wall=400,
          N=20,
          m_working=0.88,
          M=25,
          F_working=6.48,
          F_liquid=6.48,
          Cp=4200,
          M_wall=25,
          m_liquid=7.83,
          pstart(displayUnit="kPa") = 270000,
          Tout_start(displayUnit="degC") = 295.15,
          Ul=2290,
          Utp=5310,
          Uv=2870,
          U=2300,
          Tin_start=393.15,
          T1_liquid=290.15,
          Tn_liquid=296.15)
          annotation (Placement(transformation(extent={{60,-38},{28,-66}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid3(m0=7.83)
          annotation (Placement(transformation(extent={{72,-92},{52,-72}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid3(T0(
              displayUnit="degC") = 290.15)
          annotation (Placement(transformation(extent={{44,-92},{24,-72}})));
        Components.Controls.sens_PT sens_turbine_out1(redeclare package medium
            = ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={86,-8})));
        Components.Controls.SuperHeat superHeat_con(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{16,-66},{-4,-46}})));
        Modelica.Blocks.Sources.Step step1_2(
          height=-0.0229,
          offset=0.2981,
          startTime=8000) annotation (Placement(transformation(extent={{-112,
                  114},{-92,134}})));
        Modelica.Blocks.Sources.Step step1_3(
          offset=0.2981,
          height=-0.0395,
          startTime=8000) annotation (Placement(transformation(extent={{-114,
                  146},{-94,166}})));
        Modelica.Blocks.Sources.Step step1_4(
          offset=0.2981,
          height=-0.0746,
          startTime=8000)
          annotation (Placement(transformation(extent={{-120,178},{-100,198}})));
        Modelica.Blocks.Sources.Step step1_5(
          offset=0.2981,
          height=-0.1284,
          startTime=8000) annotation (Placement(transformation(extent={{-110,
                  218},{-90,238}})));
        Modelica.Blocks.Sources.Step stepT1_2(
          offset=808.15,
          height=-16,
          startTime=8000)
          annotation (Placement(transformation(extent={{-8,116},{12,136}})));
        Modelica.Blocks.Sources.Step stepT1_3(
          offset=808.15,
          height=-37,
          startTime=8000)
          annotation (Placement(transformation(extent={{-12,150},{8,170}})));
        Modelica.Blocks.Sources.Step stepT1_4(
          offset=808.15,
          height=-61,
          startTime=8000)
          annotation (Placement(transformation(extent={{-8,186},{12,206}})));
        Modelica.Blocks.Sources.Step stepT1_5(
          offset=808.15,
          height=-115,
          startTime=8000)
          annotation (Placement(transformation(extent={{-10,222},{10,242}})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT
          annotation (Placement(transformation(extent={{80,70},{60,90}})));
        Components.source_sink_ports.OutFlow_liquid_useM outFlow_liquid_useM
          annotation (Placement(transformation(extent={{-58,90},{-38,110}})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT1
          annotation (Placement(transformation(extent={{-114,26},{-94,46}})));
        Modelica.Blocks.Sources.Step stepTw1_2(
          startTime=8000,
          height=-0.5,
          offset=366.15)
          annotation (Placement(transformation(extent={{-192,10},{-172,30}})));
        Components.Controls.SuperHeat superHeat_eva_out(redeclare package
            medium = ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{38,34},{58,54}})));
        Components.main_components.Hx1counter Evaporater1(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Cp_wall=400,
          M_wall=13,
          m_working=0.88,
          m_liquid=0.298,
          M=5,
          Utp=4000,
          U=2560,
          Cp=1400,
          pstart(displayUnit="kPa") = 2000000,
          Uv=3682,
          N=10,
          Ul=3181,
          Tout_start(displayUnit="degC") = 391.15,
          F_working=2.5,
          F_liquid=2.5,
          Tin_start=368.15,
          T1_liquid=423.15,
          Tn_liquid=382.15)
          annotation (Placement(transformation(extent={{-46,46},{-14,74}})));
        Modelica.Blocks.Sources.Step stepTw1_3(
          startTime=8000,
          offset=366.15,
          height=-0.9)
          annotation (Placement(transformation(extent={{-220,40},{-200,60}})));
        Modelica.Blocks.Sources.Step stepTw1_4(
          startTime=8000,
          offset=366.15,
          height=-1.5) annotation (Placement(transformation(extent={{-216,116},
                  {-196,136}})));
        Modelica.Blocks.Sources.Step stepTw1_5(
          startTime=8000,
          offset=366.15,
          height=-2.4)
          annotation (Placement(transformation(extent={{-196,80},{-176,100}})));
        Modelica.Blocks.Sources.Step step1_1(
          offset=0.2981,
          startTime=8000,
          height=-0.1284) annotation (Placement(transformation(extent={{-210,
                  216},{-190,236}})));
        Modelica.Blocks.Sources.Step step(height=0.1284, startTime=8300)
          annotation (Placement(transformation(extent={{-210,182},{-190,202}})));
        Modelica.Blocks.Math.Add add
          annotation (Placement(transformation(extent={{-168,194},{-148,214}})));
        Modelica.Blocks.Sources.Step stepT1_1(
          offset=808.15,
          startTime=8000,
          height=-115)
          annotation (Placement(transformation(extent={{76,176},{96,196}})));
        Modelica.Blocks.Sources.Step stepT1_6(
          offset=0,
          height=115,
          startTime=8300)
          annotation (Placement(transformation(extent={{74,132},{94,152}})));
        Modelica.Blocks.Math.Add add1
          annotation (Placement(transformation(extent={{108,156},{128,176}})));
        Modelica.Blocks.Sources.Constant const3(k=366.15)
          annotation (Placement(transformation(extent={{-130,32},{-120,42}})));
        Modelica.Blocks.Sources.Step stepTw1_1(
          startTime=8000,
          offset=50,
          height=-2) annotation (Placement(transformation(extent={{-138,-52},{
                  -118,-32}})));
        Modelica.Blocks.Sources.Step stepTw1_6(
          startTime=8000,
          height=-10,
          offset=40)
          annotation (Placement(transformation(extent={{94,24},{114,44}})));
      equation
        connect(preHeater.OutFlow_liquid,outFlow_liquid1. A) annotation (Line(
              points={{-77.84,28},{-80,28},{-80,14},{-92,14}},   color={28,108,
                200}));
        connect(preHeater.OutFlow_working, sens_pre_out.inFlow) annotation (
            Line(points={{-69.48,47.8},{-69.74,47.8},{-69.74,56.2},{-72,56.2}},
              color={0,0,255}));
        connect(sens_pump_out.outFlow, preHeater.InFlow_working) annotation (
            Line(points={{-71.14,18},{-69.26,18},{-69.26,28}}, color={0,0,255}));
        connect(pump_type3_1.InFlow, tank_rou.OutFlow) annotation (Line(points=
                {{-75.28,-30.04},{-75.28,-69.02},{-44.64,-69.02},{-44.64,-68.5}},
              color={0,0,255}));
        connect(sens_pump_out.inFlow, pump_type3_1.OutFlow) annotation (Line(
              points={{-71.14,6},{-76,6},{-76,-15.34},{-87.12,-15.34}}, color={
                0,0,255}));
        connect(Condenser.OutFlow_liquid, outFlow_liquid3.A) annotation (Line(
              points={{60,-58.16},{68,-58.16},{68,-82},{72,-82}}, color={28,108,
                200}));
        connect(Condenser.InFlow_liquid, inFlow_liquid3.A) annotation (Line(
              points={{28.32,-58.44},{28.32,-70.22},{24.6,-70.22},{24.6,-82.2}},
              color={28,108,200}));
        connect(turbine_simple_control.OutFlow, sens_turbine_out1.inFlow)
          annotation (Line(points={{67.9,-0.22},{77.3,-0.22},{77.3,2},{86.2,2}},
              color={0,0,255}));
        connect(Condenser.InFlow_working, sens_turbine_out1.outFlow)
          annotation (Line(points={{60,-47.24},{74,-47.24},{74,-18},{86.2,-18}},
              color={0,0,255}));
        connect(superHeat_con.inFlow, Condenser.OutFlow_working) annotation (
            Line(points={{16,-55.8},{24,-55.8},{24,-47.52},{28.32,-47.52}},
              color={0,0,255}));
        connect(tank_rou.InFlow, superHeat_con.outFlow) annotation (Line(points=
               {{-26.4,-55},{-26,-55},{-26,-55.8},{-4,-55.8}}, color={0,0,255}));
        connect(inFlow_liquid_useT.A, Evaporater.InFlow_liquid) annotation (
            Line(points={{60.6,79.8},{35.3,79.8},{35.3,66.44},{23.68,66.44}},
              color={28,108,200}));
        connect(inFlow_liquid_useT1.A, preHeater.InFlow_liquid) annotation (
            Line(points={{-94.6,35.8},{-95.7,35.8},{-95.7,47.8},{-78.06,47.8}},
              color={28,108,200}));
        connect(Evaporater.OutFlow_working, superHeat_eva_out.inFlow)
          annotation (Line(points={{23.68,55.52},{27.84,55.52},{27.84,44.2},{38,
                44.2}}, color={0,0,255}));
        connect(superHeat_eva_out.outFlow, turbine_simple_control.InFlow)
          annotation (Line(points={{58,44.2},{58,44.2},{58,18.14},{53.5,18.14}},
              color={0,0,255}));
        connect(Evaporater1.OutFlow_working, Evaporater.InFlow_working)
          annotation (Line(points={{-14.32,55.52},{-11.16,55.52},{-11.16,55.24},
                {-8,55.24}}, color={0,0,255}));
        connect(Evaporater1.InFlow_liquid, Evaporater.OutFlow_liquid)
          annotation (Line(points={{-14.32,66.44},{-10.16,66.44},{-10.16,66.16},
                {-8,66.16}}, color={28,108,200}));
        connect(sens_pre_out.outFlow, Evaporater1.InFlow_working) annotation (
            Line(points={{-52,56.2},{-48,56.2},{-48,55.24},{-46,55.24}}, color=
                {0,0,255}));
        connect(outFlow_liquid_useM.A, Evaporater1.OutFlow_liquid) annotation (
            Line(points={{-58,100},{-74,100},{-74,66.16},{-46,66.16}}, color={
                28,108,200}));
        connect(step1_1.y,add. u1) annotation (Line(points={{-189,226},{-178.5,
                226},{-178.5,210},{-170,210}}, color={0,0,127}));
        connect(step.y,add. u2) annotation (Line(points={{-189,192},{-180,192},
                {-180,198},{-170,198}}, color={0,0,127}));
        connect(stepT1_1.y,add1. u1) annotation (Line(points={{97,186},{104,186},
                {104,172},{106,172}},
                                    color={0,0,127}));
        connect(stepT1_6.y,add1. u2) annotation (Line(points={{95,142},{104,142},
                {104,160},{106,160}},
                                    color={0,0,127}));
        connect(const3.y, inFlow_liquid_useT1.u) annotation (Line(points={{
                -119.5,37},{-119.5,49.5},{-105.2,49.5},{-105.2,40}}, color={0,0,
                127}));
        connect(step1_5.y, outFlow_liquid_useM.u) annotation (Line(points={{-89,
                228},{-68,228},{-68,104.2},{-48.2,104.2}}, color={0,0,127}));
        connect(stepT1_5.y, inFlow_liquid_useT.u) annotation (Line(points={{11,
                232},{42,232},{42,84},{71.2,84}}, color={0,0,127}));
        connect(const1.y, turbine_simple_control.u) annotation (Line(points={{
                80.5,45},{80.5,28.5},{67.45,28.5},{67.45,12.19}}, color={0,0,
                127}));
        connect(const2.y, pump_type3_1.u) annotation (Line(points={{-105.5,-3},
                {-105.5,-13.5},{-86.8,-13.5},{-86.8,-22.48}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-54,18},{42,-26}},
                lineColor={28,108,200},
                textString="ORC余热回收系统")}));
      end WHR_ORC_system;

      model WHR_OSORC_system
        Components.Controls.sens_PT sens_pre_out(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{-46,54},{-26,74}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(m0=2.1)
          annotation (Placement(transformation(extent={{-90,14},{-110,34}})));
        Components.main_components.Hx1counter preHeater(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Uv=800,
          Cp_wall=400,
          N=20,
          Cp=4.2e3,
          m_liquid=2.1,
          Utp=1300,
          M=50,
          M_wall=10,
          F_working=2.66,
          F_liquid=2.66,
          m_working=0.6,
          Ul=1200,
          U=1000,
          Tin_start=312.15,
          Tout_start=365.15,
          T1_liquid=367.15,
          Tn_liquid=358.45,
          pstart(displayUnit="kPa") = 2800000)
                          annotation (Placement(transformation(
              extent={{-10,-11},{10,11}},
              rotation=90,
              origin={-63,48})));
        Components.Controls.sens_PT sens_pump_out(redeclare package medium =
              ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
              extent={{6,7},{-6,-7}},
              rotation=-90,
              origin={-61,22})));
        Modelica.Blocks.Sources.Constant const2(k=50)
          annotation (Placement(transformation(extent={{-100,-14},{-90,-4}})));
        Components.main_components.Tank_rou tank_rou(
          L0=0.65,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          Aera=1,
          p0(displayUnit="kPa") = 277000,
          T0=294.15)
                    annotation (Placement(transformation(extent={{-38,-22},{-18,
                  -2}})));
        Components.main_components.pump_type3 pump_type3_1(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          n_pump=0.7,
          p0(displayUnit="kPa") = 277000,
          pstart_out(displayUnit="kPa") = 277000,
          T0=310.15,
          v_s=0.8992e-5)
          annotation (Placement(transformation(extent={{-56,-36},{-88,6}})));
        Components.main_components.Hx1counter Evaporater(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Cp_wall=400,
          N=20,
          M_wall=13,
          F_working=3.01,
          F_liquid=3.01,
          Cp=2200,
          M=25,
          m_working=0.6,
          m_liquid=0.272,
          Utp=4000,
          Uv=3662,
          Ul=3181,
          U=2500,
          pstart(displayUnit="kPa") = 2800000,
          Tin_start(displayUnit="degC") = 353.15,
          Tout_start(displayUnit="degC") = 444.15,
          T1_liquid(displayUnit="degC") = 593.15,
          Tn_liquid(displayUnit="degC") = 359.15)
          annotation (Placement(transformation(extent={{-8,66},{24,94}})));
        Components.main_components.Turbine_simple_control
          turbine_simple_control(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          pout_start(displayUnit="kPa") = 270000,
          V_s=1.2e-4,
          pin_start(displayUnit="kPa") = 2800000,
          Tin_start=453.15)
          annotation (Placement(transformation(extent={{74,10},{94,30}})));
        Components.Controls.sens_PT sens_turbine_out1(redeclare package medium
            = ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={108,-6})));
        Modelica.Blocks.Sources.Constant const1(k=40)
          annotation (Placement(transformation(extent={{92,42},{102,52}})));
        Components.main_components.Hx1counter_liquid Gas_oil_Hx(
          Cp1=2.34e3,
          Cp_wall=400,
          V_liquid1=0.7,
          V_liquid2=0.7,
          N=20,
          F_liquid1=2.66,
          F_liquid2=2.66,
          M1=20,
          M2=3,
          M_wall=3,
          T2_liquid(displayUnit="degC") = 808.15,
          U1=1100,
          U2=1237,
          Cp2=1.409e3,
          T1_liquid(displayUnit="degC") = 359.15,
          T1n_liquid(displayUnit="degC") = 593.15,
          T2n_liquid(displayUnit="degC") = 413.15)
          annotation (Placement(transformation(extent={{-12,112},{24,138}})));
        Components.main_components.Tank_liquid tank_liquid(
          rou=800,
          Area=0.05,
          L0=0.5,
          Cp=2.2e3,
          T0(displayUnit="K") = 643)
          annotation (Placement(transformation(extent={{92,94},{112,114}})));
        Components.main_components.Pump_liquid pump_liquid(m0=0.272)
          annotation (Placement(transformation(extent={{88,72},{68,92}})));
        Modelica.Blocks.Sources.Constant const(k=0.272)
          annotation (Placement(transformation(extent={{64,98},{74,108}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid5(T0(
              displayUnit="degC") = 290.15)
          annotation (Placement(transformation(extent={{26,-44},{6,-24}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid5(m0=5.31)
          annotation (Placement(transformation(extent={{86,-46},{66,-26}})));
        Components.main_components.Hx1counter Condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Cp_wall=400,
          N=20,
          M=25,
          Cp=4200,
          M_wall=25,
          Tout_start(displayUnit="degC") = 295.15,
          Ul=2290,
          F_working=4.71,
          F_liquid=4.71,
          pstart(displayUnit="kPa") = 277000,
          m_working=0.6,
          m_liquid=5.31,
          Utp=5310,
          Uv=2870,
          U=2300,
          Tin_start=402.15,
          T1_liquid=290.15,
          Tn_liquid=298.15)
          annotation (Placement(transformation(extent={{74,2},{42,-26}})));
        Components.Controls.SuperHeat superHeat_con_out(redeclare package
            medium = ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{16,-18},{-4,2}})));
        Components.Controls.SuperHeat superHeat_eva_out(redeclare package
            medium = ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{48,52},{68,72}})));
        Modelica.Blocks.Sources.Step step1_2(
          height=-0.0229,
          offset=0.2981,
          startTime=8000) annotation (Placement(transformation(extent={{-124,
                  192},{-104,212}})));
        Modelica.Blocks.Sources.Step step1_3(
          offset=0.2981,
          height=-0.0395,
          startTime=8000) annotation (Placement(transformation(extent={{-126,
                  224},{-106,244}})));
        Modelica.Blocks.Sources.Step step1_4(
          offset=0.2981,
          height=-0.0746,
          startTime=8000)
          annotation (Placement(transformation(extent={{-134,260},{-114,280}})));
        Modelica.Blocks.Sources.Step step1_5(
          offset=0.2981,
          height=-0.1284,
          startTime=8000) annotation (Placement(transformation(extent={{-122,
                  296},{-102,316}})));
        Modelica.Blocks.Sources.Step stepT1_2(
          offset=808.15,
          height=-16,
          startTime=8000)
          annotation (Placement(transformation(extent={{-20,194},{0,214}})));
        Modelica.Blocks.Sources.Step stepT1_3(
          offset=808.15,
          height=-37,
          startTime=8000)
          annotation (Placement(transformation(extent={{-24,228},{-4,248}})));
        Modelica.Blocks.Sources.Step stepT1_4(
          offset=808.15,
          height=-61,
          startTime=8000)
          annotation (Placement(transformation(extent={{-20,264},{0,284}})));
        Modelica.Blocks.Sources.Step stepT1_5(
          offset=808.15,
          height=-115,
          startTime=8000)
          annotation (Placement(transformation(extent={{-22,300},{-2,320}})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT
          annotation (Placement(transformation(extent={{68,148},{48,168}})));
        Components.source_sink_ports.OutFlow_liquid_useM outFlow_liquid_useM
          annotation (Placement(transformation(extent={{-48,142},{-28,162}})));
        Modelica.Blocks.Sources.Step stepTw1_2(
          startTime=8000,
          height=-0.5,
          offset=366.15)
          annotation (Placement(transformation(extent={{-182,20},{-162,40}})));
        Modelica.Blocks.Sources.Step stepTw1_3(
          startTime=8000,
          offset=366.15,
          height=-0.9) annotation (Placement(transformation(extent={{-158,120},
                  {-138,140}})));
        Modelica.Blocks.Sources.Step stepTw1_4(
          startTime=8000,
          offset=366.15,
          height=-1.5) annotation (Placement(transformation(extent={{-162,164},
                  {-142,184}})));
        Modelica.Blocks.Sources.Step stepTw1_5(
          startTime=8000,
          offset=366.15,
          height=-2.4)
          annotation (Placement(transformation(extent={{-186,90},{-166,110}})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT1
          annotation (Placement(transformation(extent={{-104,64},{-84,84}})));
        Modelica.Blocks.Sources.Step step1_1(
          offset=0.2981,
          height=-0.1284,
          startTime=8000) annotation (Placement(transformation(extent={{-200,
                  226},{-180,246}})));
        Modelica.Blocks.Sources.Step step(height=0.1284, startTime=8300)
          annotation (Placement(transformation(extent={{-200,192},{-180,212}})));
        Modelica.Blocks.Math.Add add
          annotation (Placement(transformation(extent={{-158,204},{-138,224}})));
        Modelica.Blocks.Sources.Step stepT1_1(
          offset=808.15,
          height=-115,
          startTime=8000)
          annotation (Placement(transformation(extent={{86,186},{106,206}})));
        Modelica.Blocks.Sources.Step stepT1_6(
          height=115,
          offset=0,
          startTime=8300)
          annotation (Placement(transformation(extent={{170,156},{190,176}})));
        Modelica.Blocks.Math.Add add1
          annotation (Placement(transformation(extent={{118,166},{138,186}})));
        Modelica.Blocks.Sources.Constant const3(k=366.15)
          annotation (Placement(transformation(extent={{-138,94},{-128,104}})));
      equation
        connect(preHeater.OutFlow_liquid,outFlow_liquid1. A) annotation (Line(
              points={{-67.84,38},{-70,38},{-70,24},{-90,24}},   color={28,108,
                200}));
        connect(preHeater.OutFlow_working, sens_pre_out.inFlow) annotation (
            Line(points={{-59.48,57.8},{-59.74,57.8},{-59.74,64.2},{-46,64.2}},
              color={0,0,255}));
        connect(sens_pump_out.outFlow, preHeater.InFlow_working) annotation (
            Line(points={{-61.14,28},{-59.26,28},{-59.26,38}}, color={0,0,255}));
        connect(pump_type3_1.InFlow, tank_rou.OutFlow) annotation (Line(points=
                {{-65.28,-20.04},{-65.28,-19.02},{-33.4,-19.02},{-33.4,-17}},
              color={0,0,255}));
        connect(sens_pump_out.inFlow, pump_type3_1.OutFlow) annotation (Line(
              points={{-61.14,16},{-66,16},{-66,-5.34},{-77.12,-5.34}}, color={
                0,0,255}));
        connect(const2.y, pump_type3_1.u) annotation (Line(points={{-89.5,-9},{
                -80.75,-9},{-80.75,-12.48},{-76.8,-12.48}}, color={0,0,127}));
        connect(sens_pre_out.outFlow, Evaporater.InFlow_working) annotation (
            Line(points={{-26,64.2},{-16,64.2},{-16,75.24},{-8,75.24}}, color={
                0,0,255}));
        connect(const1.y, turbine_simple_control.u) annotation (Line(points={{
                102.5,47},{102.5,33.5},{88.3,33.5},{88.3,20.7}}, color={0,0,127}));
        connect(turbine_simple_control.OutFlow, sens_turbine_out1.inFlow)
          annotation (Line(points={{88.6,13.4},{98.3,13.4},{98.3,4},{108.2,4}},
              color={0,0,255}));
        connect(Gas_oil_Hx.outflow1, tank_liquid.inFlow) annotation (Line(
              points={{23.64,119.8},{112.82,119.8},{112.82,108.2},{108.2,108.2}},
              color={28,108,200}));
        connect(tank_liquid.outFlow, pump_liquid.inFlow) annotation (Line(
              points={{97.2,98.4},{97.2,99.2},{88,99.2},{88,82}}, color={28,108,
                200}));
        connect(pump_liquid.outFlow, Evaporater.InFlow_liquid) annotation (Line(
              points={{68,82},{48,82},{48,86.44},{23.68,86.44}}, color={28,108,
                200}));
        connect(const.y, pump_liquid.u) annotation (Line(points={{74.5,103},{
                79.25,103},{79.25,85.9},{77.9,85.9}}, color={0,0,127}));
        connect(Gas_oil_Hx.inflow1, Evaporater.OutFlow_liquid) annotation (Line(
              points={{-12.36,119.8},{-12.36,100.9},{-8,100.9},{-8,86.16}},
              color={28,108,200}));
        connect(sens_turbine_out1.outFlow, Condenser.InFlow_working)
          annotation (Line(points={{108.2,-16},{90,-16},{90,-7.24},{74,-7.24}},
              color={0,0,255}));
        connect(inFlow_liquid5.A, Condenser.InFlow_liquid) annotation (Line(
              points={{6.6,-34.2},{6.6,-20.1},{42.32,-20.1},{42.32,-18.44}},
              color={28,108,200}));
        connect(Condenser.OutFlow_liquid, outFlow_liquid5.A) annotation (Line(
              points={{74,-18.16},{86,-18.16},{86,-36}}, color={28,108,200}));
        connect(tank_rou.InFlow, superHeat_con_out.outFlow) annotation (Line(
              points={{-22,-8},{-4,-8},{-4,-7.8}}, color={0,0,255}));
        connect(superHeat_con_out.inFlow, Condenser.OutFlow_working)
          annotation (Line(points={{16,-7.8},{30,-7.8},{30,-7.52},{42.32,-7.52}},
              color={0,0,255}));
        connect(Evaporater.OutFlow_working, superHeat_eva_out.inFlow)
          annotation (Line(points={{23.68,75.52},{35.84,75.52},{35.84,62.2},{48,
                62.2}}, color={0,0,255}));
        connect(superHeat_eva_out.outFlow, turbine_simple_control.InFlow)
          annotation (Line(points={{68,62.2},{74,62.2},{74,24.2},{79,24.2}},
              color={0,0,255}));
        connect(outFlow_liquid_useM.A, Gas_oil_Hx.outflow2) annotation (Line(
              points={{-48,152},{-30,152},{-30,131.5},{-12,131.5}}, color={28,
                108,200}));
        connect(inFlow_liquid_useT.A, Gas_oil_Hx.inflow2) annotation (Line(
              points={{48.6,157.8},{48.6,144.9},{23.64,144.9},{23.64,131.24}},
              color={28,108,200}));
        connect(inFlow_liquid_useT1.A, preHeater.InFlow_liquid) annotation (
            Line(points={{-84.6,73.8},{-67.3,73.8},{-67.3,57.8},{-68.06,57.8}},
              color={28,108,200}));
        connect(step1_1.y,add. u1) annotation (Line(points={{-179,236},{-168.5,
                236},{-168.5,220},{-160,220}}, color={0,0,127}));
        connect(step.y,add. u2) annotation (Line(points={{-179,202},{-170,202},
                {-170,208},{-160,208}}, color={0,0,127}));
        connect(add.y, outFlow_liquid_useM.u) annotation (Line(points={{-137,
                214},{-39.5,214},{-39.5,156.2},{-38.2,156.2}}, color={0,0,127}));
        connect(stepT1_1.y,add1. u1) annotation (Line(points={{107,196},{114,
                196},{114,182},{116,182}},
                                    color={0,0,127}));
        connect(stepT1_6.y,add1. u2) annotation (Line(points={{191,166},{114,
                166},{114,170},{116,170}},
                                    color={0,0,127}));
        connect(add1.y, inFlow_liquid_useT.u) annotation (Line(points={{139,176},
                {59.2,176},{59.2,162}},color={0,0,127}));
        connect(const3.y, inFlow_liquid_useT1.u) annotation (Line(points={{
                -127.5,99},{-111.75,99},{-111.75,78},{-95.2,78}}, color={0,0,
                127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-40,38},{56,-6}},
                lineColor={28,108,200},
                textString="OS/ORC余热回收系统")}));
      end WHR_OSORC_system;

      model WHR_ORC_control
        Components.Controls.sens_PT sens_pre_out(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{-70,30},{-50,50}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(m0=2.1)
          annotation (Placement(transformation(extent={{-90,-12},{-110,8}})));
        Components.main_components.Hx1counter preHeater(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Uv=800,
          Cp_wall=400,
          N=20,
          Cp=4.2e3,
          m_liquid=2.1,
          F_working=4.3,
          F_liquid=4.3,
          Utp=1300,
          M=50,
          M_wall=10,
          m_working=0.88,
          pstart(displayUnit="kPa") = 2428700,
          Ul=3200,
          U=2500,
          Tin_start=296.45,
          Tout_start=365.35,
          T1_liquid=366.15,
          Tn_liquid=358.45)
                          annotation (Placement(transformation(
              extent={{-10,-11},{10,11}},
              rotation=90,
              origin={-71,22})));
        Components.Controls.sens_PT sens_pump_out(redeclare package medium =
              ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
              extent={{6,7},{-6,-7}},
              rotation=-90,
              origin={-69,-4})));
        Modelica.Blocks.Sources.Constant const2(k=50)
          annotation (Placement(transformation(extent={{-108,-26},{-98,-16}})));
        Components.main_components.Tank_rou tank_rou(
          L0=0.65,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0(displayUnit="kPa") = 270020,
          Aera=1,
          T0=296.15)
                    annotation (Placement(transformation(extent={{-50,-92},{-18,
                  -62}})));
        Components.main_components.pump_type3 pump_type3_1(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          p0(displayUnit="kPa") = 270420,
          pstart_out(displayUnit="kPa") = 270420,
          n_pump=0.7,
          T0=310.15,
          v_s=1.3e-5)
          annotation (Placement(transformation(extent={{-64,-62},{-96,-20}})));
        Components.main_components.Hx1counter Evaporater(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Cp_wall=400,
          N=20,
          M_wall=13,
          m_working=0.88,
          m_liquid=0.298,
          M=5,
          Utp=4000,
          U=2560,
          Cp=1400,
          pstart(displayUnit="kPa") = 2000000,
          Uv=3682,
          Ul=5181,
          F_working=0.6,
          F_liquid=0.6,
          Tin_start=391.15,
          Tout_start(displayUnit="degC") = 433,
          T1_liquid=808.15,
          Tn_liquid=423.15)
          annotation (Placement(transformation(extent={{-6,30},{26,58}})));
        Components.main_components.Turbine_simple_control
          turbine_simple_control(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          pin_start(displayUnit="kPa") = 2828000,
          pout_start(displayUnit="kPa") = 270000,
          V_s=1.9e-4,
          Tin_start=453.15)
          annotation (Placement(transformation(extent={{58,-30},{88,4}})));
        Modelica.Blocks.Sources.Constant const1(k=40)
          annotation (Placement(transformation(extent={{140,-36},{130,-26}})));
        Components.main_components.Hx1counter Condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Cp_wall=400,
          N=20,
          m_working=0.88,
          M=25,
          F_working=6.48,
          F_liquid=6.48,
          Cp=4200,
          M_wall=25,
          m_liquid=7.83,
          pstart(displayUnit="kPa") = 270000,
          Tout_start(displayUnit="degC") = 295.15,
          Ul=2290,
          Utp=5310,
          Uv=2870,
          U=2300,
          Tin_start=393.15,
          T1_liquid=290.15,
          Tn_liquid=296.15)
          annotation (Placement(transformation(extent={{62,-54},{30,-82}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid3(m0=7.83)
          annotation (Placement(transformation(extent={{74,-108},{54,-88}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid3(T0(
              displayUnit="degC") = 290.15)
          annotation (Placement(transformation(extent={{46,-108},{26,-88}})));
        Components.Controls.sens_PT sens_turbine_out1(redeclare package medium
            = ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={90,-56})));
        Components.Controls.SuperHeat superHeat_con(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{18,-82},{-2,-62}})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT
          annotation (Placement(transformation(extent={{68,40},{48,60}})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT1
          annotation (Placement(transformation(extent={{-112,10},{-92,30}})));
        Components.Controls.SuperHeat superHeat_eva_out(redeclare package
            medium = ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{34,22},{54,42}})));
        Components.main_components.Hx1counter Evaporater1(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Cp_wall=400,
          M_wall=13,
          m_working=0.88,
          m_liquid=0.298,
          M=5,
          Utp=4000,
          U=2560,
          Cp=1400,
          pstart(displayUnit="kPa") = 2000000,
          Uv=3682,
          N=10,
          Ul=3181,
          Tout_start(displayUnit="degC") = 391.15,
          F_working=2.5,
          F_liquid=2.5,
          Tin_start=368.15,
          T1_liquid=423.15,
          Tn_liquid=382.15)
          annotation (Placement(transformation(extent={{-44,30},{-12,58}})));
        Modelica.Blocks.Sources.Constant const3(k=366.15)
          annotation (Placement(transformation(extent={{-128,16},{-118,26}})));
        Modelica.Blocks.Sources.Step step1_1(
          offset=0.2981,
          startTime=8000,
          height=-0.1284) annotation (Placement(transformation(extent={{-148,98},
                  {-128,118}})));
        Modelica.Blocks.Sources.Step step(height=0.1284, startTime=8300)
          annotation (Placement(transformation(extent={{-154,62},{-134,82}})));
        Modelica.Blocks.Math.Add add
          annotation (Placement(transformation(extent={{-108,92},{-88,112}})));
        Components.source_sink_ports.OutFlow_liquid_useM outFlow_liquid_useM
          annotation (Placement(transformation(extent={{-54,54},{-34,74}})));
        Modelica.Blocks.Sources.Step stepT1_1(
          offset=808.15,
          startTime=8000,
          height=-115)
          annotation (Placement(transformation(extent={{144,80},{124,100}})));
        Modelica.Blocks.Sources.Step stepT1_6(
          offset=0,
          height=115,
          startTime=8300)
          annotation (Placement(transformation(extent={{184,62},{164,82}})));
        Modelica.Blocks.Math.Add add1
          annotation (Placement(transformation(extent={{92,62},{72,82}})));
        Components.Controls.sens_PT sens_turbine_out2(redeclare package medium
            = ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={66,16})));
        Modelica.Blocks.Continuous.PI PI1(
          initType=Modelica.Blocks.Types.Init.NoInit,
          T=8000,
          k=0.65) annotation (Placement(transformation(extent={{-154,-74},{-144,
                  -64}})));
        Modelica.Blocks.Math.Gain gain(k=1/7)
          annotation (Placement(transformation(extent={{-166,-40},{-160,-34}})));
        Modelica.Blocks.Math.Add add2(k2=-1)
          annotation (Placement(transformation(extent={{-186,-52},{-174,-40}})));
        Modelica.Blocks.Sources.Constant const4(k=35)
          annotation (Placement(transformation(extent={{-210,-56},{-202,-48}})));
        Modelica.Blocks.Sources.Constant const5(k=50)
          annotation (Placement(transformation(extent={{-208,-24},{-200,-16}})));
        Modelica.Blocks.Math.Add add3(k2=+1)
          annotation (Placement(transformation(extent={{-124,-66},{-112,-54}})));
        Modelica.Blocks.Math.Add add4(
                                     k2=-1)
          annotation (Placement(transformation(extent={{112,4},{124,16}})));
        Modelica.Blocks.Continuous.PI PI(initType=Modelica.Blocks.Types.Init.NoInit,
            T=1000)
          annotation (Placement(transformation(extent={{132,-2},{142,8}})));
        Components.Controls.control_unit1 control_unit1_1(
          Y0=40,
          P0=24,
          t=6000)
          annotation (Placement(transformation(extent={{136,16},{150,28}})));
        Modelica.Blocks.Math.Add add5(k2=+1)
          annotation (Placement(transformation(extent={{164,-8},{176,4}})));
        Modelica.Blocks.Sources.Constant const6(k=24)
          annotation (Placement(transformation(extent={{86,-4},{96,6}})));
        Modelica.Blocks.Sources.Step step1_2(
          startTime=8000,
          height=-10,
          offset=40)      annotation (Placement(transformation(extent={{-22,2},
                  {-2,22}})));
        Modelica.Blocks.Sources.Step step1(startTime=8300, height=10)
          annotation (Placement(transformation(extent={{-22,-32},{-2,-12}})));
        Modelica.Blocks.Math.Add add6
          annotation (Placement(transformation(extent={{20,-20},{40,0}})));
        Modelica.Blocks.Math.Add add7(k2=+1)
          annotation (Placement(transformation(extent={{-140,-40},{-128,-28}})));
        Modelica.Blocks.Sources.TimeTable timeTable(table=[0.0,808.15; 8000,
              808.15; 8001,693.15; 8300,693.15; 8301,808.15; 9300,808.15; 9301,
              771.15; 11300,771.15; 11301,808.15; 11400,808.15; 11401,792.15;
              14000,792.15])
          annotation (Placement(transformation(extent={{18,96},{38,116}})));
        Modelica.Blocks.Sources.TimeTable timeTable1(table=[0.0,0.2981; 8000,
              0.2981; 8001,0.1697; 8300,0.1697; 8301,0.2981; 9300,0.2981; 9301,
              0.2586; 11300,0.2586; 11301,0.2981; 11400,0.2981; 11401,0.2752;
              14000,0.2752])
          annotation (Placement(transformation(extent={{-74,106},{-54,126}})));
      equation
        connect(preHeater.OutFlow_liquid,outFlow_liquid1. A) annotation (Line(
              points={{-75.84,12},{-78,12},{-78,-2},{-90,-2}},   color={28,108,
                200}));
        connect(preHeater.OutFlow_working, sens_pre_out.inFlow) annotation (
            Line(points={{-67.48,31.8},{-67.74,31.8},{-67.74,40.2},{-70,40.2}},
              color={0,0,255}));
        connect(sens_pump_out.outFlow, preHeater.InFlow_working) annotation (
            Line(points={{-69.14,2},{-67.26,2},{-67.26,12}}, color={0,0,255}));
        connect(pump_type3_1.InFlow, tank_rou.OutFlow) annotation (Line(points=
                {{-73.28,-46.04},{-73.28,-85.02},{-42.64,-85.02},{-42.64,-84.5}},
              color={0,0,255}));
        connect(sens_pump_out.inFlow, pump_type3_1.OutFlow) annotation (Line(
              points={{-69.14,-10},{-74,-10},{-74,-31.34},{-85.12,-31.34}},
              color={0,0,255}));
        connect(Condenser.OutFlow_liquid, outFlow_liquid3.A) annotation (Line(
              points={{62,-74.16},{70,-74.16},{70,-98},{74,-98}}, color={28,108,
                200}));
        connect(Condenser.InFlow_liquid, inFlow_liquid3.A) annotation (Line(
              points={{30.32,-74.44},{30.32,-86.22},{26.6,-86.22},{26.6,-98.2}},
              color={28,108,200}));
        connect(turbine_simple_control.OutFlow, sens_turbine_out1.inFlow)
          annotation (Line(points={{79.9,-24.22},{79.3,-24.22},{79.3,-46},{90.2,
                -46}}, color={0,0,255}));
        connect(Condenser.InFlow_working, sens_turbine_out1.outFlow)
          annotation (Line(points={{62,-63.24},{76,-63.24},{76,-66},{90.2,-66}},
              color={0,0,255}));
        connect(superHeat_con.inFlow, Condenser.OutFlow_working) annotation (
            Line(points={{18,-71.8},{26,-71.8},{26,-63.52},{30.32,-63.52}},
              color={0,0,255}));
        connect(tank_rou.InFlow, superHeat_con.outFlow) annotation (Line(points=
               {{-24.4,-71},{-24,-71},{-24,-71.8},{-2,-71.8}}, color={0,0,255}));
        connect(inFlow_liquid_useT.A, Evaporater.InFlow_liquid) annotation (
            Line(points={{48.6,49.8},{37.3,49.8},{37.3,50.44},{25.68,50.44}},
              color={28,108,200}));
        connect(inFlow_liquid_useT1.A, preHeater.InFlow_liquid) annotation (
            Line(points={{-92.6,19.8},{-93.7,19.8},{-93.7,31.8},{-76.06,31.8}},
              color={28,108,200}));
        connect(Evaporater.OutFlow_working, superHeat_eva_out.inFlow)
          annotation (Line(points={{25.68,39.52},{29.84,39.52},{29.84,32.2},{34,
                32.2}}, color={0,0,255}));
        connect(Evaporater1.OutFlow_working, Evaporater.InFlow_working)
          annotation (Line(points={{-12.32,39.52},{-9.16,39.52},{-9.16,39.24},{
                -6,39.24}}, color={0,0,255}));
        connect(Evaporater1.InFlow_liquid, Evaporater.OutFlow_liquid)
          annotation (Line(points={{-12.32,50.44},{-8.16,50.44},{-8.16,50.16},{
                -6,50.16}}, color={28,108,200}));
        connect(sens_pre_out.outFlow, Evaporater1.InFlow_working) annotation (
            Line(points={{-50,40.2},{-46,40.2},{-46,39.24},{-44,39.24}}, color=
                {0,0,255}));
        connect(outFlow_liquid_useM.A, Evaporater1.OutFlow_liquid) annotation (
            Line(points={{-54,64},{-58,64},{-58,50.16},{-44,50.16}}, color={28,
                108,200}));
        connect(const3.y, inFlow_liquid_useT1.u) annotation (Line(points={{
                -117.5,21},{-117.5,33.5},{-103.2,33.5},{-103.2,24}}, color={0,0,
                127}));
        connect(step1_1.y,add. u1) annotation (Line(points={{-127,108},{-114.5,
                108},{-110,108}},              color={0,0,127}));
        connect(step.y,add. u2) annotation (Line(points={{-133,72},{-116,72},{
                -116,96},{-110,96}},    color={0,0,127}));
        connect(stepT1_1.y,add1. u1) annotation (Line(points={{123,90},{109,90},
                {109,78},{94,78}},  color={0,0,127}));
        connect(stepT1_6.y,add1. u2) annotation (Line(points={{163,72},{163,64},
                {94,64},{94,66}},   color={0,0,127}));
        connect(superHeat_eva_out.outFlow, sens_turbine_out2.inFlow)
          annotation (Line(points={{54,32.2},{60,32.2},{60,26},{66.2,26}},
              color={0,0,255}));
        connect(sens_turbine_out2.outFlow, turbine_simple_control.InFlow)
          annotation (Line(points={{66.2,6},{65.5,6},{65.5,-5.86}}, color={0,0,
                255}));
        connect(const4.y, add2.u2) annotation (Line(points={{-201.6,-52},{-194,
                -52},{-194,-49.6},{-187.2,-49.6}}, color={0,0,127}));
        connect(add2.y, gain.u) annotation (Line(points={{-173.4,-46},{-170,-46},
                {-170,-37},{-166.6,-37}}, color={0,0,127}));
        connect(gain.y, PI1.u) annotation (Line(points={{-159.7,-37},{-156.85,
                -37},{-156.85,-69},{-155,-69}}, color={0,0,127}));
        connect(sens_turbine_out2.y_P, add4.u1) annotation (Line(points={{71.6,
                20.8},{92.8,20.8},{92.8,13.6},{110.8,13.6}}, color={0,0,127}));
        connect(const6.y, add4.u2) annotation (Line(points={{96.5,1},{103.25,1},
                {103.25,6.4},{110.8,6.4}}, color={0,0,127}));
        connect(add4.y, PI.u) annotation (Line(points={{124.6,10},{128,10},{128,
                3},{131,3}}, color={0,0,127}));
        connect(sens_turbine_out2.y_P, control_unit1_1.u) annotation (Line(
              points={{71.6,20.8},{103.8,20.8},{103.8,22},{137.26,22}}, color={
                0,0,127}));
        connect(control_unit1_1.y, add5.u1) annotation (Line(points={{147.62,22},
                {154,22},{154,1.6},{162.8,1.6}}, color={0,0,127}));
        connect(PI.y, add5.u2) annotation (Line(points={{142.5,3},{153.25,3},{
                153.25,-5.6},{162.8,-5.6}}, color={0,0,127}));
        connect(superHeat_eva_out.y, add2.u1) annotation (Line(points={{44.2,
                35.2},{-217.9,35.2},{-217.9,-42.4},{-187.2,-42.4}}, color={0,0,
                127}));
        connect(step1_2.y, add6.u1) annotation (Line(points={{-1,12},{9.5,12},{
                9.5,-4},{18,-4}}, color={0,0,127}));
        connect(step1.y, add6.u2) annotation (Line(points={{-1,-22},{8,-22},{8,
                -16},{18,-16}}, color={0,0,127}));
        connect(gain.y, add7.u2) annotation (Line(points={{-159.7,-37},{-149.85,
                -37},{-149.85,-37.6},{-141.2,-37.6}}, color={0,0,127}));
        connect(const5.y, add7.u1) annotation (Line(points={{-199.6,-20},{
                -199.6,-14},{-141.2,-14},{-141.2,-30.4}}, color={0,0,127}));
        connect(add7.y, add3.u1) annotation (Line(points={{-127.4,-34},{-127.4,
                -45},{-125.2,-45},{-125.2,-56.4}}, color={0,0,127}));
        connect(PI1.y, add3.u2) annotation (Line(points={{-143.5,-69},{-134.75,
                -69},{-134.75,-63.6},{-125.2,-63.6}}, color={0,0,127}));
        connect(timeTable1.y, outFlow_liquid_useM.u) annotation (Line(points={{
                -53,116},{-48,116},{-48,68.2},{-44.2,68.2}}, color={0,0,127}));
        connect(timeTable.y, inFlow_liquid_useT.u) annotation (Line(points={{39,
                106},{50,106},{50,54},{59.2,54}}, color={0,0,127}));
        connect(add3.y, pump_type3_1.u) annotation (Line(points={{-111.4,-60},{
                -98,-60},{-98,-38.48},{-84.8,-38.48}}, color={0,0,127}));
        connect(add5.y, turbine_simple_control.u) annotation (Line(points={{
                176.6,-2},{184,-2},{184,-11.81},{79.45,-11.81}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-40,-24},{56,-68}},
                lineColor={28,108,200},
                textString="ORC余热回收系统")}));
      end WHR_ORC_control;

      model WHR_OSORC_control
        Components.Controls.sens_PT sens_pre_out(redeclare package medium =
              ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{-66,-6},{-46,14}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid1(m0=2.1)
          annotation (Placement(transformation(extent={{-110,-46},{-130,-26}})));
        Components.main_components.Hx1counter preHeater(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Uv=800,
          Cp_wall=400,
          N=20,
          Cp=4.2e3,
          m_liquid=2.1,
          Utp=1300,
          M=50,
          M_wall=10,
          F_working=2.66,
          F_liquid=2.66,
          m_working=0.6,
          Ul=1200,
          U=1000,
          Tin_start=312.15,
          Tout_start=365.15,
          T1_liquid=367.15,
          Tn_liquid=358.45,
          pstart(displayUnit="kPa") = 2800000)
                          annotation (Placement(transformation(
              extent={{-10,-11},{10,11}},
              rotation=90,
              origin={-83,-12})));
        Components.Controls.sens_PT sens_pump_out(redeclare package medium =
              ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
              extent={{6,7},{-6,-7}},
              rotation=-90,
              origin={-81,-38})));
        Modelica.Blocks.Sources.Constant const2(k=50)
          annotation (Placement(transformation(extent={{-218,-116},{-208,-106}})));
        Components.main_components.Tank_rou tank_rou(
          L0=0.65,
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          Aera=1,
          p0(displayUnit="kPa") = 277000,
          T0=294.15)
                    annotation (Placement(transformation(extent={{-58,-82},{-38,
                  -62}})));
        Components.main_components.pump_type3 pump_type3_1(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          n_pump=0.7,
          p0(displayUnit="kPa") = 277000,
          pstart_out(displayUnit="kPa") = 277000,
          T0=310.15,
          v_s=0.8992e-5)
          annotation (Placement(transformation(extent={{-76,-96},{-108,-54}})));
        Components.main_components.Hx1counter Evaporater(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Cp_wall=400,
          N=20,
          M_wall=13,
          F_working=3.01,
          F_liquid=3.01,
          Cp=2200,
          M=25,
          m_working=0.6,
          m_liquid=0.272,
          Utp=4000,
          Uv=3662,
          Ul=3181,
          U=2500,
          pstart(displayUnit="kPa") = 2800000,
          Tin_start(displayUnit="degC") = 353.15,
          Tout_start(displayUnit="degC") = 444.15,
          T1_liquid(displayUnit="degC") = 593.15,
          Tn_liquid(displayUnit="degC") = 359.15)
          annotation (Placement(transformation(extent={{-28,6},{4,34}})));
        Components.main_components.Turbine_simple_control
          turbine_simple_control(
          redeclare package medium = ThermoCycle.Media.R245fa_CP,
          pout_start(displayUnit="kPa") = 270000,
          V_s=1.2e-4,
          pin_start(displayUnit="kPa") = 2800000,
          Tin_start=453.15)
          annotation (Placement(transformation(extent={{52,-42},{72,-22}})));
        Components.Controls.sens_PT sens_turbine_out1(redeclare package medium
            = ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={88,-66})));
        Modelica.Blocks.Sources.Constant const1(k=40)
          annotation (Placement(transformation(extent={{80,-24},{90,-14}})));
        Components.main_components.Hx1counter_liquid Gas_oil_Hx(
          Cp1=2.34e3,
          Cp_wall=400,
          V_liquid1=0.7,
          V_liquid2=0.7,
          N=20,
          F_liquid1=2.66,
          F_liquid2=2.66,
          M1=20,
          M2=3,
          M_wall=3,
          T2_liquid(displayUnit="degC") = 808.15,
          U1=1100,
          U2=1237,
          Cp2=1.409e3,
          T1_liquid(displayUnit="degC") = 359.15,
          T1n_liquid(displayUnit="degC") = 593.15,
          T2n_liquid(displayUnit="degC") = 413.15)
          annotation (Placement(transformation(extent={{-32,52},{4,78}})));
        Components.main_components.Tank_liquid tank_liquid(
          rou=800,
          Area=0.05,
          L0=0.5,
          Cp=2.2e3,
          T0(displayUnit="K") = 643)
          annotation (Placement(transformation(extent={{72,34},{92,54}})));
        Components.main_components.Pump_liquid pump_liquid(m0=0.272)
          annotation (Placement(transformation(extent={{68,12},{48,32}})));
        Modelica.Blocks.Sources.Constant const(k=0.272)
          annotation (Placement(transformation(extent={{44,38},{54,48}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid5(T0(
              displayUnit="degC") = 290.15)
          annotation (Placement(transformation(extent={{6,-104},{-14,-84}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid5(m0=5.31)
          annotation (Placement(transformation(extent={{66,-106},{46,-86}})));
        Components.main_components.Hx1counter Condenser(
          redeclare package Medium = ThermoCycle.Media.R245fa_CP,
          Cp_wall=400,
          N=20,
          M=25,
          Cp=4200,
          M_wall=25,
          Tout_start(displayUnit="degC") = 295.15,
          Ul=2290,
          F_working=4.71,
          F_liquid=4.71,
          pstart(displayUnit="kPa") = 277000,
          m_working=0.6,
          m_liquid=5.31,
          Utp=5310,
          Uv=2870,
          U=2300,
          Tin_start=402.15,
          T1_liquid=290.15,
          Tn_liquid=298.15)
          annotation (Placement(transformation(extent={{54,-58},{22,-86}})));
        Components.Controls.SuperHeat superHeat_con_out(redeclare package
            medium = ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{-4,-78},{-24,-58}})));
        Components.Controls.SuperHeat superHeat_eva_out(redeclare package
            medium = ThermoCycle.Media.R245fa_CP)
          annotation (Placement(transformation(extent={{26,0},{46,20}})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT
          annotation (Placement(transformation(extent={{48,88},{28,108}})));
        Components.source_sink_ports.OutFlow_liquid_useM outFlow_liquid_useM
          annotation (Placement(transformation(extent={{-68,80},{-48,100}})));
        Components.source_sink_ports.InFlow_liquid_useT inFlow_liquid_useT1
          annotation (Placement(transformation(extent={{-122,22},{-102,42}})));
        Modelica.Blocks.Sources.Constant const3(k=366.15)
          annotation (Placement(transformation(extent={{-156,18},{-146,28}})));
        Modelica.Blocks.Sources.TimeTable timeTable1(table=[0.0,0.2981; 8000,
              0.2981; 8001,0.1697; 8300,0.1697; 8301,0.2981; 9300,0.2981; 9301,
              0.2586; 11300,0.2586; 11301,0.2981; 11400,0.2981; 11401,0.2752;
              14000,0.2752])
          annotation (Placement(transformation(extent={{-104,110},{-84,130}})));
        Modelica.Blocks.Sources.TimeTable timeTable(table=[0.0,808.15; 8000,
              808.15; 8001,693.15; 8300,693.15; 8301,808.15; 9300,808.15; 9301,
              771.15; 11300,771.15; 11301,808.15; 11400,808.15; 11401,792.15;
              14000,792.15])
          annotation (Placement(transformation(extent={{-8,114},{12,134}})));
        Components.Controls.sens_PT sens_turbine_out2(redeclare package medium
            = ThermoCycle.Media.R245fa_CP) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={58,-4})));
        Modelica.Blocks.Math.Add add4(
                                     k2=-1)
          annotation (Placement(transformation(extent={{120,2},{132,14}})));
        Modelica.Blocks.Continuous.PI PI(initType=Modelica.Blocks.Types.Init.NoInit, T=800)
          annotation (Placement(transformation(extent={{146,-4},{156,6}})));
        Components.Controls.control_unit1 control_unit1_1(
          Y0=40,
          t=6000,
          P0=28)
          annotation (Placement(transformation(extent={{134,-18},{148,-6}})));
        Modelica.Blocks.Math.Add add5(k2=+1)
          annotation (Placement(transformation(extent={{174,0},{186,-12}})));
        Modelica.Blocks.Sources.Constant const6(k=28)
          annotation (Placement(transformation(extent={{128,30},{118,40}})));
        Modelica.Blocks.Continuous.PI PI1(
          initType=Modelica.Blocks.Types.Init.NoInit,
          T=2000,
          k=0.5)  annotation (Placement(transformation(extent={{-182,-66},{-172,
                  -56}})));
        Modelica.Blocks.Math.Gain gain(k=1/2)
          annotation (Placement(transformation(extent={{-194,-32},{-188,-26}})));
        Modelica.Blocks.Math.Add add2(k2=-1)
          annotation (Placement(transformation(extent={{-214,-52},{-202,-40}})));
        Modelica.Blocks.Sources.Constant const4(k=50)
          annotation (Placement(transformation(extent={{-242,-62},{-234,-54}})));
        Modelica.Blocks.Sources.Constant const5(k=50)
          annotation (Placement(transformation(extent={{-236,-16},{-228,-8}})));
        Modelica.Blocks.Math.Add add3(k2=+1)
          annotation (Placement(transformation(extent={{-152,-58},{-140,-46}})));
        Modelica.Blocks.Math.Add add7(k2=+1)
          annotation (Placement(transformation(extent={{-168,-32},{-156,-20}})));
        Modelica.Blocks.Sources.Step step1_5(
          offset=0.2981,
          height=-0.1284,
          startTime=8000) annotation (Placement(transformation(extent={{-146,78},
                  {-126,98}})));
        Modelica.Blocks.Sources.Step stepT1_5(
          offset=808.15,
          height=-115,
          startTime=8000)
          annotation (Placement(transformation(extent={{28,136},{48,156}})));
        Modelica.Blocks.Sources.Step step1_1(
          startTime=8000,
          height=-25.75,
          offset=50)      annotation (Placement(transformation(extent={{-158,
                  -102},{-138,-82}})));
        Modelica.Blocks.Sources.Step step1_2(
          startTime=8000,
          height=-19.7,
          offset=40)      annotation (Placement(transformation(extent={{20,-52},
                  {40,-32}})));
      equation
        connect(preHeater.OutFlow_liquid,outFlow_liquid1. A) annotation (Line(
              points={{-87.84,-22},{-90,-22},{-90,-36},{-110,-36}},
                                                                 color={28,108,
                200}));
        connect(preHeater.OutFlow_working, sens_pre_out.inFlow) annotation (
            Line(points={{-79.48,-2.2},{-79.74,-2.2},{-79.74,4.2},{-66,4.2}},
              color={0,0,255}));
        connect(sens_pump_out.outFlow, preHeater.InFlow_working) annotation (
            Line(points={{-81.14,-32},{-79.26,-32},{-79.26,-22}}, color={0,0,
                255}));
        connect(pump_type3_1.InFlow, tank_rou.OutFlow) annotation (Line(points=
                {{-85.28,-80.04},{-85.28,-79.02},{-53.4,-79.02},{-53.4,-77}},
              color={0,0,255}));
        connect(sens_pump_out.inFlow, pump_type3_1.OutFlow) annotation (Line(
              points={{-81.14,-44},{-86,-44},{-86,-65.34},{-97.12,-65.34}},
              color={0,0,255}));
        connect(sens_pre_out.outFlow, Evaporater.InFlow_working) annotation (
            Line(points={{-46,4.2},{-36,4.2},{-36,15.24},{-28,15.24}}, color={0,
                0,255}));
        connect(turbine_simple_control.OutFlow, sens_turbine_out1.inFlow)
          annotation (Line(points={{66.6,-38.6},{78.3,-38.6},{78.3,-56},{88.2,
                -56}}, color={0,0,255}));
        connect(Gas_oil_Hx.outflow1, tank_liquid.inFlow) annotation (Line(
              points={{3.64,59.8},{92.82,59.8},{92.82,48.2},{88.2,48.2}}, color=
               {28,108,200}));
        connect(tank_liquid.outFlow, pump_liquid.inFlow) annotation (Line(
              points={{77.2,38.4},{77.2,39.2},{68,39.2},{68,22}}, color={28,108,
                200}));
        connect(pump_liquid.outFlow, Evaporater.InFlow_liquid) annotation (Line(
              points={{48,22},{28,22},{28,26.44},{3.68,26.44}}, color={28,108,
                200}));
        connect(const.y, pump_liquid.u) annotation (Line(points={{54.5,43},{
                59.25,43},{59.25,25.9},{57.9,25.9}}, color={0,0,127}));
        connect(Gas_oil_Hx.inflow1, Evaporater.OutFlow_liquid) annotation (Line(
              points={{-32.36,59.8},{-32.36,40.9},{-28,40.9},{-28,26.16}},
              color={28,108,200}));
        connect(sens_turbine_out1.outFlow, Condenser.InFlow_working)
          annotation (Line(points={{88.2,-76},{70,-76},{70,-67.24},{54,-67.24}},
              color={0,0,255}));
        connect(inFlow_liquid5.A, Condenser.InFlow_liquid) annotation (Line(
              points={{-13.4,-94.2},{-13.4,-80.1},{22.32,-80.1},{22.32,-78.44}},
              color={28,108,200}));
        connect(Condenser.OutFlow_liquid, outFlow_liquid5.A) annotation (Line(
              points={{54,-78.16},{66,-78.16},{66,-96}}, color={28,108,200}));
        connect(tank_rou.InFlow, superHeat_con_out.outFlow) annotation (Line(
              points={{-42,-68},{-24,-68},{-24,-67.8}}, color={0,0,255}));
        connect(superHeat_con_out.inFlow, Condenser.OutFlow_working)
          annotation (Line(points={{-4,-67.8},{10,-67.8},{10,-67.52},{22.32,
                -67.52}}, color={0,0,255}));
        connect(Evaporater.OutFlow_working, superHeat_eva_out.inFlow)
          annotation (Line(points={{3.68,15.52},{15.84,15.52},{15.84,10.2},{26,
                10.2}}, color={0,0,255}));
        connect(outFlow_liquid_useM.A, Gas_oil_Hx.outflow2) annotation (Line(
              points={{-68,90},{-76,90},{-76,71.5},{-32,71.5}}, color={28,108,
                200}));
        connect(inFlow_liquid_useT.A, Gas_oil_Hx.inflow2) annotation (Line(
              points={{28.6,97.8},{28.6,84.9},{3.64,84.9},{3.64,71.24}}, color=
                {28,108,200}));
        connect(inFlow_liquid_useT1.A, preHeater.InFlow_liquid) annotation (
            Line(points={{-102.6,31.8},{-87.3,31.8},{-87.3,-2.2},{-88.06,-2.2}},
              color={28,108,200}));
        connect(const3.y, inFlow_liquid_useT1.u) annotation (Line(points={{
                -145.5,23},{-131.75,23},{-131.75,36},{-113.2,36}}, color={0,0,
                127}));
        connect(superHeat_eva_out.outFlow, sens_turbine_out2.inFlow)
          annotation (Line(points={{46,10.2},{58,10.2},{58,6},{58.2,6}}, color=
                {0,0,255}));
        connect(sens_turbine_out2.y_P, add4.u1) annotation (Line(points={{63.6,
                0.8},{104.8,0.8},{104.8,11.6},{118.8,11.6}}, color={0,0,127}));
        connect(const6.y, add4.u2) annotation (Line(points={{117.5,35},{109.25,
                35},{109.25,4.4},{118.8,4.4}}, color={0,0,127}));
        connect(add4.y, PI.u) annotation (Line(points={{132.6,8},{140,8},{140,1},
                {145,1}}, color={0,0,127}));
        connect(sens_turbine_out2.y_P, control_unit1_1.u) annotation (Line(
              points={{63.6,0.8},{97.8,0.8},{97.8,-12},{135.26,-12}}, color={0,
                0,127}));
        connect(control_unit1_1.y, add5.u1) annotation (Line(points={{145.62,
                -12},{164,-12},{164,-9.6},{172.8,-9.6}}, color={0,0,127}));
        connect(PI.y, add5.u2) annotation (Line(points={{156.5,1},{165.25,1},{
                165.25,-2.4},{172.8,-2.4}}, color={0,0,127}));
        connect(sens_turbine_out2.outFlow, turbine_simple_control.InFlow)
          annotation (Line(points={{58.2,-14},{57,-14},{57,-27.8}}, color={0,0,
                255}));
        connect(const4.y, add2.u2) annotation (Line(points={{-233.6,-58},{-222,
                -58},{-222,-49.6},{-215.2,-49.6}}, color={0,0,127}));
        connect(add2.y, gain.u) annotation (Line(points={{-201.4,-46},{-198,-46},
                {-198,-29},{-194.6,-29}}, color={0,0,127}));
        connect(gain.y, PI1.u) annotation (Line(points={{-187.7,-29},{-184.85,
                -29},{-184.85,-61},{-183,-61}}, color={0,0,127}));
        connect(superHeat_eva_out.y, add2.u1) annotation (Line(points={{36.2,
                13.2},{-251.9,13.2},{-251.9,-42.4},{-215.2,-42.4}}, color={0,0,
                127}));
        connect(gain.y, add7.u2) annotation (Line(points={{-187.7,-29},{-177.85,
                -29},{-177.85,-29.6},{-169.2,-29.6}}, color={0,0,127}));
        connect(const5.y, add7.u1) annotation (Line(points={{-227.6,-12},{
                -227.6,-6},{-169.2,-6},{-169.2,-22.4}}, color={0,0,127}));
        connect(add7.y, add3.u1) annotation (Line(points={{-155.4,-26},{-155.4,
                -37},{-153.2,-37},{-153.2,-48.4}}, color={0,0,127}));
        connect(PI1.y, add3.u2) annotation (Line(points={{-171.5,-61},{-162.75,
                -61},{-162.75,-55.6},{-153.2,-55.6}}, color={0,0,127}));
        connect(step1_5.y, outFlow_liquid_useM.u) annotation (Line(points={{
                -125,88},{-90.5,88},{-90.5,94.2},{-58.2,94.2}}, color={0,0,127}));
        connect(stepT1_5.y, inFlow_liquid_useT.u) annotation (Line(points={{49,
                146},{49,125},{39.2,125},{39.2,102}}, color={0,0,127}));
        connect(add3.y, pump_type3_1.u) annotation (Line(points={{-139.4,-52},{
                -120,-52},{-120,-72.48},{-96.8,-72.48}}, color={0,0,127}));
        connect(add5.y, turbine_simple_control.u) annotation (Line(points={{
                186.6,-6},{194.3,-6},{194.3,-31.3},{66.3,-31.3}}, color={0,0,
                127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Text(
                extent={{-60,-2},{36,-46}},
                lineColor={28,108,200},
                textString="OS/ORC余热回收系统")}));
      end WHR_OSORC_control;
    end WHR_system;

    package Solar_system
      model Solar_ORC
        Components.main_components.SolarCollector_New solarCollector_New
          annotation (Placement(transformation(extent={{-46,-8},{36,48}})));
        Components.source_sink_ports.InFlow_liquid inFlow_liquid2(
                                                                 T0(displayUnit=
               "K") = 527)
          annotation (Placement(transformation(extent={{-88,8},{-68,28}})));
        Components.source_sink_ports.OutFlow_liquid outFlow_liquid2(
                                                                   m0=15.91)
          annotation (Placement(transformation(extent={{54,8},{74,28}})));
        Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(
          tableOnFile=true,
          columns={2},
          tableName="nana",
          fileName="nana.mat")
          annotation (Placement(transformation(extent={{-106,62},{-92,76}})));
        Modelica.Blocks.Sources.Constant const2(k=25)
          annotation (Placement(transformation(extent={{-4,84},{16,104}})));
      equation
        connect(solarCollector_New.outFlow, outFlow_liquid2.A) annotation (Line(
              points={{36,20.56},{46,20.56},{46,18},{54,18}}, color={28,108,200}));
        connect(inFlow_liquid2.A, solarCollector_New.inFlow) annotation (Line(
              points={{-68.6,17.8},{-56.3,17.8},{-56.3,20},{-46,20}}, color={28,
                108,200}));
        connect(combiTimeTable.y[1], solarCollector_New.I_in) annotation (Line(
              points={{-91.3,69},{-50.65,69},{-50.65,36.24},{-9.92,36.24}},
              color={0,0,127}));
        connect(const2.y, solarCollector_New.Tam_in) annotation (Line(points={{
                17,94},{24,94},{24,36.24},{13.86,36.24}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end Solar_ORC;
    end Solar_system;
    annotation (Icon(graphics={
          Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={28,108,200},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-60,80},{-60,-80},{62,0},{-60,80}},
            lineColor={28,108,200},
            fillColor={170,213,255},
            fillPattern=FillPattern.Forward),
          Ellipse(
            extent={{-30,20},{8,-14}},
            lineColor={28,108,200},
            fillColor={255,255,0},
            fillPattern=FillPattern.Forward)}));
  end Example;

  package Media2Type
    model R245fa_CP
      extends ExternalMedia.Media.CoolPropMedium(mediumName="R245fa");
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end R245fa_CP;

    model R123_CP
      extends ExternalMedia.Media.CoolPropMedium(mediumName="R245fa");
    end R123_CP;
  end Media2Type;

  package Cache

    model MB_model_inner
    //medium//
       replaceable package medium = ExternalMedia.Media.CoolPropMedium (mediumName="R245fa");
    //variables//
      Modelica.SIunits.Pressure P;
      Modelica.SIunits.Temperature Tr1(start=350);
      Modelica.SIunits.Temperature Tr2(start=353);
      Modelica.SIunits.Temperature Tr3(start=358);
      Modelica.SIunits.Temperature Tw1;
      Modelica.SIunits.Temperature Tw2;
      Modelica.SIunits.Temperature Tw3;
      Modelica.SIunits.CoefficientOfHeatTransfer U1;
      Modelica.SIunits.CoefficientOfHeatTransfer U2;
      Modelica.SIunits.CoefficientOfHeatTransfer U3;
      Modelica.SIunits.MassFlowRate m12(start=0.1);
      Modelica.SIunits.MassFlowRate m23(start=0.1);
      Modelica.SIunits.MassFlowRate min;
      Modelica.SIunits.MassFlowRate mout;
      Modelica.SIunits.Density rou_cool;
      Modelica.SIunits.Density rou_heat;
      Modelica.SIunits.Density rou_l;
      Modelica.SIunits.Density rou_g;
      Modelica.SIunits.Length l1;
      Modelica.SIunits.Length l2;
      Modelica.SIunits.Length l3;
      Modelica.SIunits.SpecificEnthalpy hin;
      Modelica.SIunits.SpecificEnthalpy hout;
      Modelica.SIunits.SpecificEnthalpy h_l;
      Modelica.SIunits.SpecificEnthalpy h_g;
      Modelica.SIunits.SpecificEnthalpy h_cool;
      Modelica.SIunits.SpecificEnthalpy h_heat;
      Real drou1dp_h;
      Real drou1dh_p;
      Real dhldp;
      Real dhgdp;
      Real gama;
      Real drouldp;
      Real drougdp;
      Real drou3dh_p;
      Real drou3dp_h;
      medium.SaturationProperties Sat1;
      medium.ThermodynamicState state_cool;
      medium.ThermodynamicState state_heat;
      //parameters//
      parameter Modelica.SIunits.Area A=4.9e-4;
      parameter Modelica.SIunits.Diameter Di=0.025;
      parameter Modelica.SIunits.Length l=10;
      Modelica.Blocks.Interfaces.RealInput u_m annotation (Placement(transformation(
              extent={{-74,-10},{-34,30}}), iconTransformation(extent={{-70,12},{-46,
                36}})));
      Modelica.Blocks.Interfaces.RealInput u_h annotation (Placement(transformation(
              extent={{-78,-36},{-38,4}}), iconTransformation(extent={{-70,-24},{-46,
                0}})));
    initial equation
      P=8e5;
      l1=0.6;
      l2=9;
      hout=4.7e5;
    equation
      //出入口参数//
      min=u_m;//0.1;
      mout=0.1;
      hin=u_h;//3e5;//medium.specificEnthalpy_pT(P,273.15+50);
      //---------------------------------------------subcooled Region---------------------------------------------//
      //-质量守恒-//
      (rou_cool-rou_l)*der(l1)+l1*(drou1dp_h+0.5*drou1dh_p*dhldp)*der(P)+0.5*l1*drou1dh_p*der(hin)=(min-m12)/A;
      Sat1=medium.setSat_p(P);
      state_cool=medium.setState_ph(P,h_cool);
      h_l=Sat1.hl;
      h_g=Sat1.hv;
      rou_l=Sat1.dl;
      rou_g=Sat1.dv;
      h_cool=0.5*(hin+h_l);
      rou_cool=medium.density_ph(P,h_cool);
      drou1dp_h=medium.density_derp_h(state_cool);
      drou1dh_p=medium.density_derh_p(state_cool);
      dhldp=Sat1.dhldp;//medium.dBubbleEnthalpy_dPressure(Sat1);
      //-能量守恒-//
      A*((rou_cool*h_cool-rou_l*h_l))*der(l1)+l1*A*(0.5*rou_cool+h_cool*drou1dh_p)*der(hin)+0.5*A*l1*(rou_cool*dhldp+(hin+h_l)*(drou1dp_h+0.5*drou1dh_p*dhldp-2))*der(P)=min*hin-m12*h_l+3.14*Di*l1*U1*(Tw1-Tr1);
      Tw1=380;
      Tr1=state_cool.T;
      U1=800;
      //----------------------------------------------two-phase Region----------------------------------------------//
      //-质量守恒-//
      (rou_l-rou_g)*der(l1)+(1-gama)*(rou_l-rou_g)*der(l2)+l2*(gama*drougdp+(1-gama)*drouldp)*der(P)=(m12-m23)/A;
      drouldp=Sat1.ddldp;//medium.dBubbleDensity_dPressure(Sat1);
      drougdp=Sat1.ddvdp;//medium.dDewDensity_dPressure(Sat1);
      gama=0.3;//需要修正//
      //-能量守恒-//
      A*l2*(gama*(h_g*drougdp+rou_g*dhgdp)+(1-gama)*(h_l*drouldp+rou_l*dhldp)-1)*der(P)+A*(rou_l*h_l-rou_g*h_g)*der(l1)+A*(1-gama)*(rou_l*h_l-rou_g*h_g)*der(l2)=m12*h_l-m23*h_g+3.14*Di*U2*l2*(Tw2-Tr2);
      dhgdp=Sat1.dhvdp;//medium.dDewEnthalpy_dPressure(Sat1);
      //dhldp=medium.dBubbleEnthalpy_dPressure(Sat1);
      U2=800;
      Tw2=380;
      Tr2=Sat1.Tsat;
      //----------------------------------------------superheated Region----------------------------------------------//
      //-质量守恒-//
      A*l3*(0.5*drou3dh_p*dhgdp+drou3dp_h)*der(P)+A*((rou_g-rou_heat)*der(l1)+(rou_g-rou_heat)*der(l2))+0.5*A*l3*drou3dh_p*der(hout)=m23-mout;
      h_heat=0.5*(h_g+hout);
      state_heat=medium.setState_ph(P,h_heat);
      drou3dh_p=medium.density_derh_p(state_heat);
      drou3dp_h=medium.density_derp_h(state_heat);
      rou_heat=state_heat.d;
      //-能量守恒-//
      A*(rou_g*h_g-0.5*rou_heat*(h_g+hout))*(der(l1)+der(l2))+A*l3*(0.5*(h_g+hout)*(0.5*drou3dh_p*dhgdp+drou3dp_h)+0.5*rou_heat*dhgdp-1)*der(P)+A*(0.5*rou_heat*l3+0.25*drou3dh_p*(h_g+hout)*l3)*der(hout)=m23*h_g-mout*hout+3.14*Di*U3*l3*(Tw3-Tr3);
      Tw3=380;
      Tr3=state_heat.T;
      U3=800;
      //-----------------------------------------------------其他方程-------------------------------------------------//
      l1+l2+l3=l;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(extent={{-62,46},{66,-32}}, lineColor={28,108,200},
              fillColor={175,175,175},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-12,24},{62,-16}},
              lineColor={28,108,200},
              fillColor={175,175,175},
              fillPattern=FillPattern.Solid,
              textString="MB_inner"),
            Text(
              extent={{-46,32},{-30,16}},
              lineColor={28,108,200},
              fillColor={175,175,175},
              fillPattern=FillPattern.Solid,
              textString="m"),
            Text(
              extent={{-46,-4},{-30,-20}},
              lineColor={28,108,200},
              fillColor={175,175,175},
              fillPattern=FillPattern.Solid,
              textString="h")}),                                     Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        uses(Modelica(version="3.2.1")),
                  Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end MB_model_inner;

    model Test_MB_inner
      MB_model_inner mB_model_inner(l=12)
        annotation (Placement(transformation(extent={{-12,-24},{68,34}})));
      h_source h_source1
        annotation (Placement(transformation(extent={{-66,-24},{-32,4}})));
      Modelica.Blocks.Sources.Constant const(k=0.1)
        annotation (Placement(transformation(extent={{-62,18},{-50,30}})));
      Modelica.Blocks.Sources.Constant const1(k=348)
        annotation (Placement(transformation(extent={{-82,-14},{-70,-2}})));
      Modelica.Blocks.Sources.Step step(
        offset=348,
        startTime=2,
        height=0)
        annotation (Placement(transformation(extent={{-88,-50},{-68,-30}})));
      Modelica.Blocks.Sources.TimeTable timeTable(table=[0,348; 1,348; 1.5,345;
            2,345; 5,345])
        annotation (Placement(transformation(extent={{-66,52},{-46,72}})));
      Modelica.Blocks.Sources.Step step1(
        startTime=2,
        height=-0.1e5,
        offset=2.8e5)
        annotation (Placement(transformation(extent={{-40,-54},{-20,-34}})));
    equation
      connect(const.y, mB_model_inner.u_m) annotation (Line(points={{-49.4,24},
              {-22,24},{-22,11.96},{4.8,11.96}}, color={0,0,127}));
      connect(const1.y, h_source1.u) annotation (Line(points={{-69.4,-8},{-56,
              -8},{-56,-8.88},{-57.5,-8.88}}, color={0,0,127}));
      connect(step1.y, mB_model_inner.u_h) annotation (Line(points={{-19,-44},{
              -8,-44},{-8,1.52},{4.8,1.52}}, color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Test_MB_inner;

    model h_source
      replaceable package medium = ExternalMedia.Media.CoolPropMedium (mediumName="R245fa");
      Modelica.SIunits.SpecificEnthalpy h;
      Modelica.SIunits.Temperature T;
      parameter Modelica.SIunits.Pressure P=8e5;
      Modelica.Blocks.Interfaces.RealInput u annotation (Placement(transformation(
              extent={{-80,-22},{-40,18}}), iconTransformation(extent={{-60,-2},{-40,
                18}})));
      Modelica.Blocks.Interfaces.RealOutput y annotation (Placement(transformation(
              extent={{60,-2},{80,18}}), iconTransformation(extent={{60,-2},{80,18}})));
    equation
      h=medium.specificEnthalpy_pT(P,T);
      u=T;
      y=h;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
              extent={{-56,34},{66,-22}},
              lineColor={28,108,200},
              fillColor={175,175,175},
              fillPattern=FillPattern.Solid), Text(
              extent={{-30,18},{34,-10}},
              lineColor={28,108,200},
              fillColor={175,175,175},
              fillPattern=FillPattern.Solid,
              textString="source")}),                                Diagram(
            coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
              extent={{-58,32},{64,-24}},
              lineColor={28,108,200},
              fillColor={175,175,175},
              fillPattern=FillPattern.Solid)}));
    end h_source;
  end Cache;
  annotation (uses(Modelica(version="3.2.1")), Icon(graphics={
        Rectangle(
          extent={{-80,60},{80,-76}},
          lineColor={28,108,200},
          fillColor={85,170,255},
          fillPattern=FillPattern.HorizontalCylinder),
        Bitmap(
          extent={{-2,-20},{84,56}},
          imageSource=
              "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"
               +
              "1VHju0Dc/sZviNssEh2IMVXjvMLNVWc8vldemV4X0r3OXQpnmxgAzhhqwgaue2ce7p62AATCA2cPx5Pnj0SaxwFBBEvY/Y7pl6YOuvHg5wc8KKSKSgabsUICawUSyEnxYkjnQkTjNm58aw7eWbAJcXsvoRgAl8bnrXnhnVCfFmkAgKkrSvHpusou09bUvOkZN3nS3StqB9w9bU3yHef0QFuFxarfuX6zrjnfKUz2jncobCqBwSIMO+2Be71cP3Bvq8kzlmLsKX11QdQBAJLcOlSemAmuP3EgXr10DPoV5fzbzk9EyPK5cGL/duBcwXcrd6C6MQS3Q0Oy29E8I2QkuaBrCrZV1GN9cQ08LgdOGdABuvbvm4IzDmfT4twfteBx6nDqCmwpEYwZSHFq0BSGmGGiRaoP7159HD666SScOrQLCIQBHQswddKJuHPcEGgqR3aKFy9cdhQ+v/UUXHREd9B+YwVBZQxH92yN968di5MGdwARwaEpcOkapIRm2NQLy6pw2WOfHpDMh7VLw2WDChoHt06+zaXz2bVxcfTs4sZpDTH7k6glJ1mC+hORW+Vsl1NRP3Gp/KqiVOexGV7HmvZZbnfHSxOGqdE987D2uY/DqiKXJ2qQGIQawgYem74cU75cAqeuYnSPIvhciUFEEjpurhJpNcHf10D+UxwwqXbtqEcgGAUROhIhFYxja2UjtpTVQVM5wnGjaRtmny1fxhqcqv61c9wpePOOM7F69WZUN0bhDxuDI5Y4Iyrk4yGTvi0OGJ+/srNuwsowtR967YuOuhnV2FIaqHZp/Psmfjq4wlVNP/CxoThgYGN9LJUEtWUAiusCuOi56Xhv3joYpo2zRvTAu9cch+N7tv71pQgBwzq1QNfCTFT7w3huxnLUBsJw6ioK0xNrKM45ijKTADBsLqvHs9OXw7JtDOnUAq0zkn6VsARAV5XEugiEumAMqR4nXLqGmGWhLhRHVrIHisLREDYQihsIRA34wwYaglHUNoYRiRnQVBWMcdQHo6gLRNAYjsMfjsO0BYCmjkUEjQOXj+mJt68+DmN6FqG6MYIpXy7FFa/OQiCWiKwVQnbufmxH54/bSw9I5q/fcD4mz95KWyobcwNR892oJZ6M2TRGEPk4eL1D5fPTXOqNLZKcR3ROc5zRJs33bn3E6rGtIfZRdcQeVBwwAABPXXISPGceCZ+uzuOMNftRMQaYQkBhHApXsGx7BcobI00+k9QiYpiFjVHjgPvOz3HAbkrHPj4VZ/dtpXy5ufpwYqRwBjRGDMzfUIoBHQpw7ojumLmuFKbYOyIoDLu8Ot/u053YuufDMIMzWd8tLHuuIeQAIkqJCzEiLjFMYaK8PhZfnXRk8lpOoLApT5IJdWrj0Z1ygh8s23HAgtA1FSrnpsdpzoRhJ1uC8hduq8HK4pno3yYbd58+FCO6t8ILlx8F47kZmLl2936Bhg5NwVnDOsPp0PH9umKs2F2DcFxAVTiKMpIABqgKR6uMhKWzIWzgu/Wl2FXtR/u8dJzYrx0enb7slwUjgs+pITPZDdO2UdYYRGGGD05dQ119FHWhGPJT3QAYagIRbK30Y9wTXyLN5wCQ2Ag+Y2hH3HHaUFQ0hHHR81+jyh8BZwymbaO0PtScZEbhDJeN7okHzhoGp67i/blr8diXy7CxvB6WJHDGpK7w7UlOZU5espMEgPUHKPeueRlIdalbKrfUrWQkWzo4m68oypcOzhblehxb1je6wrnJcb2kMXxy0I5fZksaREQuDhoc88fnfDh9KcaN7Y+C5GQ4VGVRQ7RxkyD0axIdMn0uHN23DUzLwlfLt8OyJRhnIEKyIdDJXxdevn71WnTt+ef8RH8NBzzER+MSC4sbU+NC9iNK5LQmIsxZW4xgJIYB7XPRKjNp75mtBCiMbx7RMTuQnuQFAPTs2RE053bUvnv9ty3TnaemObWTPZrykEPh8zngF1K2MKQcGzbE7UFT3GFJ6qIyVLo5/+iyd36S3z444YAFsemFy5GtmY3juuXcnOpQj3Fr/GanpqwwbWn+sKkcFzz/DT6Yvx756T48dv5IdGmR0TyzEBHyU73o2SqxUfzjpjI0REyETQsAQ2qTtY8re9+HDQuljSEsb3J9Gt2zNbwOrdmXsRmMIRAz8NlPWzBt/gZsLK1H29xUMMZQ40/MOAu3VGD60i34ZPEWmEKizB/B2pJ6rC2px4aS2kQuRQCmLbCpvAHrS+qxtqQOmyv8iFoCAANJwhFdC3DHaYdB1xQ8PX0ZJr4+B2tKakFA2KUqs3xO9bIsj+OYV8af+kTPth2N7x77c+H1+2L25Evw8fIKI9WhPJnq1E7IS/GOC5eXPlvfIFe4HJwyHA1Hb6uLvVNv4lVDyMOJyKFwttOh8lp8uwShSGJiSnEoOLpthp8xtrnZrEmE7gWZ6JifidW7qjF7XUnztggBnKQYMXFMa3XitJ8OuB774oBnqkhirZRDYPnNKh5nWLi1HD+sL8bY/u0xoE0utpQ3IBGpSdA51r763RaBb+7cp+8kKrvp+asbbnvmk5nn9mgx8+QPl/v8Edk5bMbGxQUbRES5ALjKqTTdyR4dkp40e1swjIPlHTn7qaswG7DG3vXB2gnDBqy95bOZb1X5QxNDBptYXBdIvfm9H5CT6sXIbq1xztDOuO2DH5sbr31uCnLTvPCHoli8pQwkActK7Oa7XQ4AiexKLmdCpzctG9KSWLCxFOOGdkb7vFTkpXqwtcq/n6mbAYiaAnd+vBAcidwX7XMS69A1OysRNizM21SKhdsqE6Pwnt/u6Ty0fw4MxhnAf34qI8HnVHDZ6J7ITPbgpW+W495PFiJs2tBVZbdPUx5ulZH0/oq14cjh/TScOLLlQZJ4E76+AxXAUgDQz3sYPTq31Gsi1H1jbfRGw5bHSIIXAKkMxU6H8q9Uh/ri8BbJu7ZcfSyq6xMuXD89fQl+Oup+mZTmXGcKCUp0QxzRrQBup475G0pRE4qBs0SUOQGI2dTvi43+VMZY7cGszgHPVEHLggXZkgjNm0RNm3j4bs1ugHGc0L8NfG4dTUE8JhTa6Ez698aFh64+FZ2GD8SmFyaGmBlc0ilVvyHfqx+Z5tWHZ3gdw7NdOPaqfoXTNYeQq16bdHAbGMD0+87ACaNaQ9PsmkGdku9LdisTVUUJlNaFMPmzJYgaBkb1aIU0r6PZba9bYRacuoaVu6qwsrgGYIDd5CKTpGvgnMOpKPA2na1ky0RayLnri1FSE0B2ig9tc1P/rR+gkIApJDwOFYWZCVEv3l4FyxIgMJj2Ho99ljjNXSReEHK/fbDmz/a8KJHmYEiHFhjVswhldQE89+1KhOM2NEUpzvJo59ZdOfLVDGZGMPd6fNbkbX4wMWNVKbqOfxBvTPuWaarad3uD/WJtzP46aonTJcGjc7YrxaHcluXWxnROoRuTddr+7r0XiqXPXInbrzy7+T6OZBXEaCMDGQmP/BScMrAjYqaJ79cXA3J/PUASy46bdnYkbh7U+hzwTGVEDWictQahecOJiOByaMhN9QAgdCnMQk6yB9uq/GAKb9DBd6mqjvgfuH/FtDtQkVhJB5teAIBb3jqocvhVrJ0yEfrEZ+z2qa4vVhrR82wuj1y1uwY7qvxom5OKLnlpmL+lApwzdCpIBwCsL6lDMGbBoe4dNPTEaepgLLFuAQAhJMAYyhoj2FHdiNY5aeiYl44ZK3buVwYGNGeaAoAMnws5qR7ETQu7a0LN8mZNI7BLUzCmexHyM5ITjumSMLxrAQAg2ePAeYd3R2MknlhTWTbmrS/G9ko/uhVmwutyYv6GYuysCYBzBh9nHz59bNGio2aspJlT/rw/5R/B+qo47pi1cUR9RL5mSVmkcA6PxrerjH2a49G/KsrybFiwo1pdUiVzFW4x95kPK6maHpk8+pS6maXb8c5to5Cs6uAMu+KMGoSUuQXpPuSkecEA5Ka4f7HZTaAkIUWLqCUOdGm4Hw6IVIsWrcXga6ZCtNFb7y1xQs25fHRPXHt8f4RjMTz2+U/YWesH4wSHgs1FqXq5JYF6AN2ufBqcseTqsHV63GJhxlDuUFm1g7NGj1sPjWyfZExfvkUWv33fPxLCPv3szsi9a0lESWYbOPiRwbiF6kAUXQsyUZiZAmwuR6rbgY55ifPrNpfXJTZV901zyNCscuyBbQuACDHTwvaKBozq0QZdWqRBVTlEsxd/wncyL8UFKYGqQAxd8jOQnexBMGqgqjEMlyMRuFgViMKwBNrlpODlS49CdqrvF3XJTPLgobOG7/fZ5M8X4+Z35iG9ab+wrD6CmCXAGYsKhm9OeX+DxGe3/+VydqUmIS5ANiGbmiLCFcbVmC2H7A4aI3YH45oluYeDqSQBg5ijUYo5zy6debFBMAGgfaYHDpWXL60IbLSIcn/aXoXnvl6Om08dhLtPTxhqvl1XDNaUfYcITkOgZdw4uGb1AyLVlW8twBFD2qkLK+tbEpqC7Ihh3MD2uO3UwSAC7p02H+8t2AxJe3LcKcuXPDvD/8Xc+3HC80BZwABnaBUwxMMkWTqAGOMswIA6hI2anXXhMgW+0uSzHt3V8YoX1sYtuzgvWQksevwm481pH+OCM4/5Sxv7tgW1wHe3Qx/3sCtuE1waR4rLgahhoi4QAYjQLjcN7fPSETdN7KwOAE3ZoZqxD0n2QNcTqcykpMQ6CkDfNrnISnKjwp8w+0oiFKR78d7VYxNm/SnTUZiVBJdTQ1l9AA3RKC4f0wNXHdMX93/8I96YuwEVjRF8umgjerfJbX52VqoHrXPSEDNtbC6pgdmUaCcQMzFv7W4AEpF4wmMjzeuEpnAIIsXp0NIDDZHmmfCvwjG9CpA//nEwRjvDDNUgViSEREjIVgBrRWyPDPdRYwmwSPYpCQZSGFgNADx50RHo3+36gO+s4cs52BFR08Jj05chy+fEhUf2wdMXH4nzpnyJJTuqmiy3BCHREpY8qHU8IFK1yPSAQCmMsQJqSkF2XK/WePT8kfC5dEz5cimen7U2YU7fW+CdGDsEJ4xMJDVUwOBUuRk1aJ7gKCRQJgg5gpBDoEQqMACMQeysCzUAqGuM8lfQ6qopLqfjT5b8j+PL1aVoO2FK29KQMZRIoE1OMlpnJ6O4Noh15XUAgJaZyUjxutAQiqGqMdKUJHMvIk3ZYW0BxO2mgw32yWZUUhuEEBL5GUnITvGgwp+w1oGAFmnJ6F2UC38khmS3A1nJHgAMgagJW0gM7pCPVtmp6JKfDs4Z6sMxTHpvAbwOranTSFx0RA9MPv9wVDWGMf75b1DeEAbngGlLhGImAAVbKxsgbBtdCjOQm+xBSUPIEYqZ47u3zv6+zYTnG/9qObfNTIEtqdoQoTtsSfkAJANTGCMdDIyIsbBhJluCtMSYRQoD8zs1Te7xtu/frR1w/L3gkLsYEiTxx0zc/uFCJPncOP2wLnhuwmhc+vIsrNxdA3AGBl7w/eXD+dkPvHfQpqsDItWy0lpwzgqFlC2JgCHtc/HoeSOQl56MF2csw4OfLkbMpqbplsAYM3VVLXO6nc3rqUW3nAe3g28++qmPzjRs4VSZmtIYi/WKCeoqiDq6dG1kKG62iFtSkYwyOZDJJFfR0YVxJxyOMw6kAr+BUXe+jO/uuwStJkzpUB017zel7Opx6DhvWDekJ7nxxZKtqA5EAcaQm+qBoihoiMRQF4klSMUYtCajRMRM5DE3pUSsyYdPU5Q9WVlQ6Q8jYpjwuBzISnY3uQQAIILHoUBVGaKGjahpNxl8EiHvlg14XYmBZXttAFJIMIUhbhNiVpMbkSCETRtNfnuoj8RRH441E58xAArDws1lWL6tAn07tMAJfdvi+e9WwxByTGkwcmOqS30GRedVXXbzGXjp0r9GM/jhsQkAYACY9nvX7tnKGHf/28rCreVgyl5V2+l2Q1PVMlhxE5Q44b4mFMcd0xYg0+fCyB6t8dDZw3DZK7Owuy4IAcp7ePF2pyQ6aGmhD4hUfkPC51DzBVEqCOhSmIWOBQmPgldmr0UwbjVPswDAgChnssqt82ZSTZv7I5LIICeDIhRk1QRDg2M26y9AHSRRS0vIdKspBzgYDAfHzqwU53yfi//6yRgHiN5XPItMp8K3Vxmtks95/IK4JU6zBDowABOP7IGLx/RCQyiGd+dvhC0IYAz5aYkN3Rp/DKFoU+phDuhaorEjMQNAQp2LxBKWJodDbZZLXTCGYNREiwwHclM82MsqINmlQeMccdOClISkJhIl0kYDjj2hHsaeYDwGQO49svNnhx0kgvn2mtn3PKmsMYJX56xBn/b5uHPcENQGo/hg0SY9IOmmqGUfkT6w67sL1+z8ZOKjH1dvrwvQt48d+N7gv8PSFcuQnexiRz8yQ2UOnyMQNT2pHj21IRJt1xA3WrrOfCQPpKQzhiyQ+nosbDTnoHNrHAxUxYlHJEhPtAXDtspGPDV9GQZ1aoFeRTkoSPdhd20QKmMFoaiZC+DAPQiacECkalJlMomggzFMX74dZ6zbjWFdW2JwhzysLdnf/E9AiKta/Z74lpNvfR7ba6scM3ZFLwmbdAIRtZZELSSRzgDiHH5J2KBxttHBlXVcwXqPwrYVpCm7qxtsHJiTzF70vvRp+Jy6ti0QbrnNH++/WdJIW9JQW6KdlJKnuB24ZHQP3HDCABARnvh8MRZvqwDjCefN7JTEhm5pTSPiTTORU1WQ1EQcfyChzklB8AcT71OcOhSmQEAiFLcQisYBJCEryb2fRcOhq1A4gxAykcqtyc/RMC1oicPNAFCz+T7hWLuPuk30i3B+UJNpmfYcjMBAjOPzZdsxrNN6nHd4dzxwzjCoCsMnS7YqcUv2t0n0CteFL9zZGJ2ncmV+u0ueWnbWkHYVuxsi9PZ14w5KOxSc/ySSnHqbkY99f6GAkgOJTIFQKkDpVZF4EgNSbcmbMq0m5OHVlRWxGJpJpXAGhbMGySkEidSmSsPt0HB07yI4dQ1vzVmLJdurAM5gCJG9uSqUa9riv4NUCBvQM5JSDWGBM4bKxgjembcew7q2xPkju+HrlTtR2hBpbl/OENSEEeZNo/D83QFwhvSQxS6ziHUGEim+vA5tg2XLWdlex2y3ky3v0yq7vsIfFt/cdwEaAZQdrNo3YXfIhMekXvVR8Y5FaEtECkmCU1PRoygHVxzZA2cO6w5bEt74biWe+WZlIqsuY9AUhgyfCwCwsbwh4UvHWCK3utMB25bYXZuIQLWkQHFtwgye7HZAVRiEzRAxbOzxQctIdu+3e8ibZpaYJWBL2Tz72UKAg0HlCQ8W0eQGluTSMGFkV7TJT4SUC5LoVpgOgJCR5MZDZw9H2LDAGEPcMPHxwi1YvKMCjDHUR0zcNu0HpPucOLpvO7xw6VHoVpCJt39Yjy3Vfs2wZS9TohdnNNGWtPDTFTvOtiT9oWN2/ggilgWTRIYlcKkt7XRizRaKPSDOEGKMhaWkVAKclhB5mHkXFi3ZiMEDOiOJEyRkWAEFZXP+fmBIhzycPaIbGsMxfLBoS8JYwxkA5nJoaoYpDp7W86dJRURgxz0I07KzsKfqnOObVbvw08YSDOnaCuOHd8X9n/2E5sMJgIYkryeico4SAC6PDwpnyY2BaDoTiXMaLEkIGXYLW9Ix8UB0AA+w+q01u8oUxkvyz31ia8tM35KSqmD51FuOxrDuXQ6KEIKWRExaHkuiJYGUVJeOUd1aYtxhnTGkUwGyU7yoC4Rx/0c/4vW56xAz97cUiURGS+yoDjTrU05VgduhwbAFSgIJVxoShN0NCVJ5XDpUzmFAIGZa2FJWh8EdCyBs+asJNg0hYYOa97n2vYAoEcgImcgkdNPJg5GV4t37fSLvHZLdDpw3ott+95WSsHh7eZNKCJT7o5jw0re4YWwDrjymN248ZTBOH9oJs9fsxrQfN2PJ9gpELaEZgrUPGOT5rZD8/xT9W+eAgN3zd9SU2GSnOxS2VGHKbMaRFrXERSDIJI3dkJHkmVvujz4Vs+WxTk1rfcb973mnfLUsDACaywlbyqgMC/8eQmqKgrOHdUaS24lnp/+En7ZXNC8diEhriJlpfyBc7g/jT5Pq22Xrga9fA427IqV5zcSA6mAMT89Yjt7t8nHakI54e94GlDSG9nTCcF0oZtlNcUkuDqgcDR6VP2tx0c4SLJ+AbEmUDVBbU1CHvXY0ibjNbH+l/+FIdfguVdUOmhCIAFuSoKZ91quO6YObTxoEl0OHPxzDvxZtxEszV2H+5nJYgvaLADZsiSe+WIJ563Zj3saSRJJQIqR4nPA5dcRNq2kdlRh1/ZE4pBRI8+hwORRETBuWkHhq+jJsrWjAVyt2/GrYvkvj0Pbx0tgXnDE4dAXgDLtqgrjvox/RNi/hykSS0L11Fo7o1hr+cAyfLN6MUNwCY0DcsPCvn7ZgX4IyxlAVjOPeTxZh5e5qXHR4dwzqWIAJY3rhpIEdce0bs/Hegg2QTDGDcWkeRE6hItgITeF+IWUlgF4KZ3MipnV7hlvvFbfYOGKUlOZxRbYt3rLN0bGgjAEIG3bRl1vrUiRRGABiVhwulVkcFBBN9e/XIQ9H9WqDivog3pi7AYYp9pIKgK4pXivwR1wR/hj+NKk27igD0TqWdPbL7n090MEZZq8rxaodlRjcqQBnD+2ER75cCgBwqkrj6onD7eeWFOOh94Ctr1wFANUPTp72YPv8DPbA4h2OoGF6PTovqAjEu5pSFtlSKSSyW4CQI4Bsp6aGIkxgcOf2B00IxABBJMBIghgykr1wOXTMW7cLz3y1DDPXlSBqNBld2M+DZRnmbSzDvI2lYIwn4gWJkOJxwOXUUe0Pwx+NN/fb2mAMcdOGz+2Ax6GhLpRI/LKuvAHryur3MyIAiXsRERyJs4VhWAlZq4oCQWjartg7c4YNC8/PXAXsSelsS1x2TG8c0a016kIx3PPRwkS+P86aZ9Wf788wBoRNC1N/3IxvV+3AqQM74OaTBqEoNx2+poBLIrIj0V9h+AHg3J6tcMOFx5res6bUWVKCMZ5Dn9zKci582k8QYUlIrQxEuiMr9X2NU60lCJzDneLhLltKNAK4eFARbj+qi510w5d+y7Th0BVcfWwfZKf68P68tdhQVr+fgYsB0Bglw/gvINW/ftyM+cu2qrYQ3v2+kIR0nwvpSS4ISSASTS40DLaQ4bz+feSsnzbhIQDdr3gGTo1pz2wo02qXbZM5Dl1+cM2JDUP7dqgDjFXvfj0POkXx8co6bWcD8wYiZqpPRSjk8+FgemtpnIEBwpSMiAi2nTA2TF2wEZ8t2Qamqb95nlSTbr7fZ13z0+HSVQQicUQNq/n7QNRAzLCQ7HYize1CMYX2EunnMxQDDFNCSglNVQDOETESZXPoaiKltpnwMt9jvk+UZ29Zif+aQ+2e5/2WVBg4BxqiFt7+YSOO6VWEotx0iL08EvJn4Y0HitGD+4MNvIXUwqwwAJAUnXtc/GSqxlnQpfHlhi2LPU59e3RXMQrad5tRGYjWORntaOFllYYNFAMY0LUdWF5v6Trj7kiiFgymmegtualeeJ06GvdJNUcghE0rFT4V20tr0bYpWvpA8KdJtbTKBCC5YPDse6BbqseBe04fjA4tMvHZ4s14Zc66xMjOCKZAGMffj/5dE2fDbq6NgjE6ypbsBoJm1pqwjnvq6xDhK8OWIm4IMohIujU14lCVIJGMx3V1ZlKKVlt3EBszx+uAwnm4NBiLS2H74mai47p1bb+c738UjDF0Lkwkfqn2R5oMA4mJoSEchz9qID89Cfmpbqwq3ms+/zUE4gaEJOiaCsYZIk1lczSRyGwK3dCVg38q0p4IWkXh8DR52u/JleHV1YZ2GW7DlhJrD9LzenbKA059HF5VLDGheE3b3rwzYIvQh2cFjrl1znllUdtukeU1v1mzGZuev3IxgMUAsG9azIEdsoBjb0dcyDAAxC0b9328CJ1bZGJY15a46qgeePSLZTDlnoPhGAhwoIZA8uBMvH+aVF63CiLiIcPWmi22BJw2sD3GDe2K7RX1mPzFUjREDPCmvSqVwTQVDSnehAnalAIK5zkSchjAYUtC0DSx/+QMxCwbMVtCIbJdCi+Nm9amg9SOiXIAUEBxIrIAINRkiUvzOv5YWrKfwaGpKGxKu1bZGGpW2QCGYMxEbSiKNrlp6NU6G18ljqr5t4iZEpYg6KoCTUnsVwGA26FCYaz577wk1y/OzDoYIAIcqoIklxOSJGJNM6Wu8PL7R3UJ72gI4/+1d9ZxVhX//3/NnLh9t3uXZWFZuktKBEnBwAJFRQxM7MDuACXEQBFFVESxAIuS7u5Ylu3O23Fi5vfH3SVMYlf9fH/7fDzW2L33zDlzzvvMzHve79f7vnpsL1JW0DI+fKk3iD2Krus6I0lpd/4cy0E4B0dRVrlubN2MkTaTwTjhuso1SderGEcwuPhJ2MJigUtfh00iTo8S2kfMKnNh5k87MOvO4Xjg0p7YeqwUy/blnSjAJ1Aq6RvnISvvynq5hnM2KrORQGcgniAoR2hB2K9lPJ68ujcY53jtu83YcqwYpFbjTwBlTSOsBZmZJ6WtDJTCKEvZ3qD2IwM1CAQmSiEK4JLOmUnh3AhOBEqISAFBEEQmSaJbrE9XDQCH2wsCBDgnAXCg0hsyKqv57MOgOAcsBhFxEaFZcUGVB6yu6iIBfIqGsmofAILWqTEht/qfLvZDRqOqKmwmGYl2E5ye0LnZTTJkkcLpC01tmsaGQRJCL6b6hcMoi7CbZagqC0VjAHArWuWIy/oHN2zbX6+tVTsVHJMCV7oD2ssMjIETBkCvPRXUBt7Xukg5owQ+s8U0waeoG04cxOtBSkp88dEKr65zLoASfL0lEx1TY/DQqN544soLcKiwCgU1HhACiFQQgzwHi1buqJdrOGejSo2wQ9UZrfHUiLrGEGE14PEreyM1NhwL1+3Hd1uPnfaW5+BcFhE8tcJEWjAbkoZVuSxpjWCPI21irEKMVaA04Bf3lTukHK9XkiUjjTJZZLPEZEogR9uM+TaDAEc93sjm0RHQOfdnV3tcOjS4vEGAMzSJsMIkCfDrDGc+BnBEWoyhTVxw5FY4T9SPAkLZt/lVIZHL5EgbTLIEd1D90+O7/Aq8ioZomwlJUVaUOEN7XmZjKEerxhtaYMdF2iEKoXKvZ3GyZ4TFIMFikhFQNbj8Sl2RvkplyAvo26P9+TdwGgL8KjeoOouq68/QP2v/rbNaRyoBFSg4QVCQiJVrpyddngoB4Fd1vL1sN/q3S0X/9ml4/IqeePjTVahX92Ut52xUkiAC0AmpTeSJCzOjfWosVE3D6gN5cPiV33mV+G82YI544xGVGJ7g8+vNxKBHPFoVlLNrIIFzY7VPtXJIoqLrxBP0VjBNOhwbEXmk2qvqh987+0Juf4UgATKBj4CXAQTlLi8CqoYWiVEItxjhd3rPfBrIORIjLAi3GuHxB5FX7jz9IWccOWUOAEB8hBURFgPcAeVPCgWERjZfUIMhSkZshA355Q74AgqsBgNsJgmOOpFNuxlmWYRf0UDq06p4SADUbpLh8gVC9Z5CHsPyBinDbpDBdaw0SOL9skBFTrhOOUtyB/W7LQbJMq5faxgNMgqq3PhhZxb8ig4OCn7qNYfSJU7vSkpQXO3BvpwydG2RhG7pCbCbZFR46l/45ZxXt55gEN6gyhnjjFCCY6VOfLl2PyRRwL0je6JFXPhp4TEhZxPlp2W2WozwKtqljLPlKtN/rvEFl1R4gt9XeJWvNIZPwMkczsiH7iD7piqg/FruqLqza6JV7PXguVed+CMyYkzYOf2uQJhJPgAKVLr9cPsVJEbbEBNmOiE2dEZwoGlsOGwmGceLq3G0sBK/rZ59vLQGQVVFjN2CxAjLb+//afiDKly1a7zkCBuOlTpR4w3CapJhNxtrR65QAuKpRafrjdoXpsUoI7e0BtUePyghukRRaqiNJKlXvnscvq8ePvjlQGGm8/P2065Kj5tlFLCHc06MkohbBnfClHEXY9IVPWtDtOrO8tQykjooKD8tVYQx9MpIwogeGfAHFcxZuQcVnjMoOHEOnLNRFTsDKHUFGOOh+a7OGN76ZTd+3p6Ftk1iMXF4V5ikUxL1OIRKTzAa0inijQKFxiEAMHHABPCQ7jO4DiAAwA3AzQGqMZ5c7VOfXZFZ0uVQmbteO+GFq/uDDHueK7p+kBDCy51+FFa4YDcbkBhuPbuDkVChAIDiSHE1Sly/kYumBIcLq5Ff5kSYxYC02LA/TaEnhMCvqrXVAIGYMDM8fgUuXxBGWUS0zYAatw+cc9hMMsLMxr8sy3OupETZIAoC9udXweELghL4CSUlUv07HIHhL6LprTPa3rBKfcI0evfHXxwqWlHlY29zEHOVJ4DHP12Dt3/ejjcXb0O1LwiAcFVj+mnyboIJvoBirhv+OedoGmPH86P7IDbcijkr9mD+hiO/mx386/lUXr8GzsE451rdCRXXeDF1yTZ0TU/ALRd3wJ6cMny89gAoIWDgKPP6o2mYGTqvACExoAJglrDOx8mdOudauFF0G0UhoDAEfaquBlQtKAsCBMLH+FV2t8Z5jCOottYY31af97FNqzSYRr8BCppHCPFXuv3mncdL0Tk9AU1j7Tg1avzvkEUBLWujGTJLqqHp+mn7RoQQFFR5cKiwAi2So9EyMeovD61qDOWOUJhThNWAoMZQ6fKhdUoM4sMsyCp1wK+oMEgiLMb6izI5ASVoFh+qGnKs1AGmcwgicRokqUwSRJx9sdK/IaDAo+ijfKr+EjsxkeW6LNCtOuNxy/flNV2xL//EhrdIcYRDzZFFfnLv0mZDvtufpPNQxws0VCTjovZNsed4Md5duis0Ta7dX1QZU4g4Fr9ueLJeLuGcjUpVNVDCOcC1uqeCCBRrjxThpa/WY9ptQ/DoqJ7YmlWEQ0U1ICCgVDBqORU4nBOazrBvnoArJB13AACq/6Addu3roIRQQsh4xrlFZUxqgJcxbEYBlJBCSohDUzXzgYIKACT00J/hG5lzjrToMHRMiwVnDEeKqv9wahfQdBwurMblAPpkJCHCbETNn9SKUnWOT1bvQ7jZiG83HoI7qKCs1si6psZi8fYsvL90J6xGCcdLq+vdrW6QBLRKjISm6zheGspVJISWWARSSUSK+pT435/rRvt7psOt6HYQeERCyiTKj8qCtMhuEJZ4gkoPn8peVTlpChAmE34ozCQ8WHqg7Hi3gU2xA8DezCx0vPlTsOSQKg/nDMM7NsWEwZ1R5fbjsc/X4WiJ47SoCp0xBSM7oUVaQr1cxzkb1ZVto9Ah2qq9tKHA6dNObprpnGPBpiMY068N+rZtiruGdMaD81ZD54BBglXdNhnzf9h0xu3Emg2wGaSccp/yvDeoJ1oMdJ9X1c5INOZsMIkiOHg5AS8CQeLRkhooqopWyZEwyRIC6hlsDHKO9k1jkBYbAYc3gGMljj92QHCOAwUV0JmODmmxaBJjQ01e5R9+llCCXw8VYdOxJQiqOnTOUVwdWke1SY1DUGN44ov1oJQgqLF6dVJwzhFrNSI9IRxVrkCo/E4o+uN490SLK9+p4Vg93oNYuxV92kYhr1L5RhTkNQRsb6zMqzpHBgKlQQtMkH/aVhHcWeynqQIhLNYo5CYZWHly93jseCMkTJNfUgO+8QXYrp9m9zAVNqOMB0Z0R1SYBZ+t3oc1h/J/N80ziMTn97rRJD6qXq7jnI1qUOeWGH9ZX/bihqne3/7NZjIi0maGxx/E+kMFoZR4AnBGbLe8PF8sdPn/sChTj1s/gxZkxGsKGAK63xbUWZRbZXFO1Z8KDgOhPFzTWZuAL1Cv0z8AiLZwtIoxur85EDgGkO45ZU5UuHxomRiJ5AgLjpU5/37OTYD2KdEQRRFFVZUoqHL98dSOEGSW1MDpi7IU4AAAXjhJREFUDSLSZkJGQjj25oZGxj/DV1dom3EU1ka6J0RYYDcbUFzjAWF/9t3zMDLOcUFGAprGhWNvTjmKHT5QQiBQdnDB9mKd//gUyPTzF9SsY9eeTdh4qBLxkdYmjgC7BxxrFIMwL7NazO3RNAIGk0hkl8ufbsShZraAG7KVLX759Eovs5buwKrthwXOmS2kmahhY2YxBnZKQ0yYGSZJgktXTusVgyDV+JX/QOpHtzbNQPo/y4U4y2nTakoJru+TgTZNovHzjuP4cXcOeO3DGND1sO+OlomEQAOAjhNegUC4WBqUEn0qTTjiL02klLYMOlkfDSyNcxbFQGzgMPHaSZhOqYYlz33ywXebcMeVveutI3bOvB87R7yq2yMMe1Rdvb64xoNjxdXo3y4Vg9ql4ljp3r91q0tUCJUCBbDlSCGqPIE/NkQC5Fd6kF/hRKdmCWgWG/G3z/6pfy6odIIxBpvFALtJQnE1/uL77ITEGeMc7Cz2ZQgh6N+2CQySjL3ZZXD4ghAEOGIt5o1eSa13MZj31uZgbJ82xkV7c8b6Nf0imZIkiYhfS6KA7KrKyGo/v92rsDGcc+Q76dIwyfE+Rjyd16FVK+ybegMAYH2eE+vznJJfJxFAyIE2Z+VeXHVBBi5s0wT92yThh13Zp3WYT2NuWOrPk3nORtWhRRIw5HmYRaHaUzs14oyjQ5No3DW8G9x+BR+v2gtvXXQ3AIlSe/NIq0GkJLANQHHQAIGShGqf9pXGWXsCbuJcDxURIuCyQDWTQCWnoiNU/YnAIHDx7Xc+oYvW7az3cg2JMUbIgrAsoOj3u4Nq0o6sElzUPg192zbBh2sO4O9yhzQGzN9wCJUeP+as3AdVZ3/44BFCUOnx4YUv12NAh1Qs35d75vtglGDj0WLMWLIFVW4/iqr/Zg+NUGw4WIC5K3ZjT245qtyBMzIGzkPKSt3TE+DzB/Ht1kwwxiAKwhGbLOwyCEC9yroC2JhVBkJIa4XxCwACgdLPclcWHGo5rGlYvtP3TpDhWh7yFkNT9E5BjV+QFht5nbOq7ET4X+uYGGgM8oHyGpuqayHHUI0HX244hJfG9scjl/fEtqxSlLn8tds8YBKBkxsFnF8dzpOc38aGRKEDJ2JbKSEY1SMDKTFhWLB2P37Zk3PaHo3OeWSh02vlnDsBwCYQUAK/m5Ig07mbEnKMEBSBkyyDRBzhJqmzovERUJggElTJlHwTLtO371hRyZKE8nq+pUBcuB0iaFapJ3gQ4EnL9+TgtsGd0aV5PBLCLCio9vz1s084ft6TG7pu0L98eBkHFu3MxuJdxwFCz/itT2oFOB+bvy6UFlKbbvKnn6cEh4qqMOHD5aEab2c6uHCGFvFhaJkUjeNlDuzNrwyF9BC6Z9+7dzkueeYzHK3n/nf6g6CU2HXADlCoDOUomINy/3M3KoyM4pwLIiX5koBjQY33Ujl6ORS9a0BnP9Udo6DaCQAWzvXwk9dC8N2WTIzt1xa9WiZhcIcm+HzdYUCgACGqSZaqtGB9mdR5yj6LodyaGgIwzkLFxm4a2B5ev4LP1h6EL6CdfsM5bDIR7GKt0OR9A9Lx4KAW1Qk2w4QIk+HCRLtxwMXpEdfE2aXvNIbuRc7g0ApPkBgFuinZKt/YMzH23kij+SAWP4Ki76bU8y0FrmiVgO1bDvkMVNhCiYCtx0uxO6cUzePC0a1ZHM4k04EQAgJ6ZjW3KQFAz2HVw8FAwM/09hESkgU5i4YIpRjasTnCLCZsO1aMCrcXhFJdoHwjueo1/PLyTfXT6aeQGG1DQqS1VBaEEg4OibDrrNc/86orqL3AODeKhFTEW+V7OiWEjRYJDnDODV6FJ/iVU5xIkgjBIIVx8JMy5ITjUHE1Zi3bBVGkuHlAB8SGW8BDw56ba3qlpR69pudlVBQEsiiWgJAAOEevFgloGhuGoK7/YWAnA7d6NDXGoYRc6g/ceDnuue4y1sSWmFlR5Tqm6zxuXbbr9WJX4Fu/pl8iUJpvMYgPpIQZRuUcKfvFaqawW1hC+/FPxYRdOrHeb+pzd46AKTUOBpGuI4R4XX4Vv+7LgySKuHVAW9hMMs5kRXLGowHwB0mPZ/ilU/51JvxWIfcvP8s5mkXbMPbC1tB0HSv25kDXGGSBZkdY5C0J4cZz7OG/ZmjbZOR9cG+mTPk8QsB9Ovp7Ff0JxnkkJeBGEbM/ujTpJ3+AGTiIhYDqZgkVRvHkY+zyBeALBGPAcVKil4RktzWdQdMYujSLR4eUaIBziKKQ3715bEHPFnH1dh3nZVQpdhFRRhTLlJSDAou3H8PCDQcRaTXhhdF90SUtBpydFpVl5Jwk6Kd4qlJvnoESX4UlItx8U4UvON+rafcZJCk6wWbaFGURJ8sSzS71BEebU6OmLz9etnBHkW9ZYcB6q7NSQFVe/UZWAECzcDNSbIZdBoFsBDhWHchDpdOD3m1T0aFJzGmVzv/Pwjn6tEpG84RI7M0pwabMYoAKMAhkSd7jF2UPbxHdIM3OuX8UYsa+xm0Sf98qCq+JhBwQKSkRKD1mlujUGLs0c8jcHF7lDVzMOVIpISUGUTxik09G6eggYCDxADlh+ZxxXN6tOZ66pg8IIXh/2S5sySoBKAVjrKja4a3Jr3DU23Wc15qqQ2oUBELyi/eXZBNCmhbUePD45+sQZ7fgoo7NMHXcxRj3zk/Ir/aE1gwcgs5YE8V9sl5sdVCBR9fbuhV9us55JBBSrfWpeiuPorzGgQjOINYpB4YEIMW2943pJ4ye9X29pnMDwMEP7gEuf63GbhK/opQO2pVTThdvy8Stg7tgZNfm2JhZb+JB/1mMkojLu7eAQCk+X3sIBVVuCILglkTxJ9PTy1ngy8fqvc3UW98AIdTkDbKuOqEVXZOsk3Od/pkCoVavovuax1oqVF3Tspe0gOH27FU2HY9Rops6RFuzA1w44TQJeHwQ7aZUzrkAhEbdnulxeO2Gi5AQacPHK/ZgyuLt8ARDSxNKSMm2mXcqs7/fhAn1VP73vEaqMEnGB5dd5ANBQWhhTpFX5cYds5dje2YhLurQFO9PGILEMHNtcC2HQEmbHZOvEV/6fAUAQBREgNBqgJxQCPUENTgDaiRnNIwS0WkQ6W6bQZprlYUnrLIwxiIL07ZlFcHvd9b7zQWAOLsMm0FYJ1KSH1Q1fLXxCLyBIC7v0QIZcX9cRvTv4LxWc+93c7aQBsWJ0jaM/WVdYa7X/vxJUe86TYtTf87qPBnH8E5NMahjUxRUOLCs1ukiCdgbbSS7U2zyWR3vTPEqBD4FaY6g/kWlT1uzqcC5tNQVeLHGq/SPsshmrz8oXZIeT+6f5kPm69cUHXjy4lnvDoiavvnD59T1028FAEz7aiX4p7cJhJBWrPZaOiZHYc5dQ9EyORrfbjqMSfPXweEPnihOLhDkY/gLmDCq/rZnzmukmvv4aJDhL3BThDG/TqGNEoLM0ho88ulqzLhlEIZ1TceLo/vi8fnrUO0NIKjp3a99Z1mczngRAFhlQON6hUTJFsppokCRRwlyKCV5MiX5kiyVRBlo6cuXdK6qoZr+xIKdnYKadm2hmwZUxuo1A7iOFrFWpEYacxcfKFuh6fT27cfLsO5QPoZ3Sccdgzvh8S/W/0Vi4e/hnCPObkas3YjMUgdOrKs5IFKCZvHhyEgIh0EUUFTtweGimtobf0JIAZQCTaJsaB4fAYlSFFa7kVnqOKE/CAACIYgOM59Ita9ru8rth0/RzuBEAbtRxm0Xd4TdYsTHK/fgeLkTIiUBmyzOO7Ii23HtTZ3rNYqiDr9OYBRpAkCMjPNoXefxBOilMH28u8JTRAmyMiuyDlBK9n6558uDZlnKT7bwasPVD6meuY8CAD5YcRhzfj0Sq2m8K+ccLeLDMfXmgWiXGo8lWzPx6GdrUOk9uaVAOBglyCMW6RzWtX/OeRkVIQS4+lUQ0OyQqDEIR8hztO5QAWb+sBVzJl6G/u1TEWk1osobgMZ4So1PaaLqrAgACjweIBBwJUZETYi2mtSOkdT/2UtzdM63gxCClWt3YNDN7+MJjUTWBNQbaoLawxrjKUTQj9UUOg5f/ern+ObJG+r1Bm+YegeOj5umJdmNc3Jq/Jc4vIGkaYu3oldGMq7p0wZzVx/EgcJKUCr8vpzoH0AJweOX98Q1vVviprd+wOrDRSAkpIFx//AuuGNoJyRE2SAQCrc/iDX78zBp/hocKXGAkJD+xLgBbfHIZT2QGhMGSkJS0V+sO4gpP+xAWW31kV6tkvHO7UMQYTWdolDLMHvZTryyaBv+zrPBOUP/tino374pjhVXYdayXVA0BqNE90SZpe8iBjfHwufH12tf1xFpNULX9G3xNtMVnmBwkC/I+jHw5gw8TmNoCqApBx9EAOZXmcsR0Aq9Af2ZKpd/Ud0xSr1+SJSmKDpLAweaRNtwQaskOH0BTFuyFbnlLpBTtTwI9wkEuSaRot6E1FEPRd8MogBBoLlc5X6Am0+eMEWLpCgIlGLRliPILneGXMgcVp2TdI/CNgMAFj4LALwYqCkGToiIEELQ5753sXVPJonr23JYrsv7uKqz3owTiRJoImEpGx7vSZ9cmV/vm8AA0MQu4dIWETsmby5ZpDF+z4bMUizddQxjLuyARy/rgXs+Wgmvop7Rpq0sELRICEdytB0t4yKw+nAhBEoxcVgnPDumHwKKhkWbj6LC7cfAdk1w+QUtEWkzYtw7PyGn3IXbBnfA6zcOgCgIWLknG56ghgHtUvHgFT0hSwIenLcGmqajVUoUOjSNhaqxEzobnIfc9n8XaM85kBJpw1NXXQCLUcZXGw4jq9wNgRIui1h4JLOielj3pJOFz+uZgtn3AID7jpkLN1klvvmXfcVGlyrEeoNqM0XX2wlE6BVQtd4a4/E65+EEXJVkuSTSLpwIxHYFNViNUgYAGwjB7twybDyYj0Gd09EsPhxrDxed1iYBrRZFqVgQxP+WUdkkGYSgwEf1Ss7QJHSDOFIirRjaMQ0BRcGv+/OhMwZCKRgH1XXSDruzMePLVXhgzMDfHfOJD5bg1QmXIu22t1q+tqXw9qDGrtc4EgBApuSgWaJvp9gtX/e9fyGb+vQErGuAm7z17YnYP2YKsxnlL4I6vyKgqElzV+/HoE7NcVXvVli0/Ri+33YMRDiD6AScdLOH5BVCmbo3DGgPgQLPLliLD1bsQ1DT0TUtFp/efyn6tmmC3hnJKHNm4to+rWExSHj5qw1488cdCGo6hnRoio/vvQRXXtASn645gO1ZJRAFAQDBxkMFeOTTVVA1BgaOvCo3/tb/zhmuuSADPTKSsC+7BF+sPwTGOUwi3RZulL43JYdh6eQ7G6CngXEPT0N0TLgYbjUzWaTs8Tuv5gD8CKmO5cWNfX310A5NZq84WNzEr+vpOmNdKQSDRTQctkmh7IafflqPEc/+AL1FVFteq5Zd7Q5g8Y4sDOnaAgPapeLzDUdRJ+RKAEgUeZEmqRLgqE91rvNOM4u3GRBvM1RIFLknfsk4LmqXgs7NE3GspAZ78ypOFnYmgMZ5+z6X9jbM33j6e+/deV8C/V7C0r0Ftugbp15f5gku9Kn6w5zzBIMQyn0RKVlTM//RDyyyXo2jn+LhG/s2yI0GgMvapaB08ugtZonOphR8zcFCfLR8NywmGXcP7YzYMPMZ5QTyP/g/QaAwyiJ0naGoyo2gqgHgKHZ44KzVpTPIImRBhM0oQ9V0ZJVUwuUPIqhoyCyphiegIjHKjq7N40/bHFM1hmKHB4dLqnAgvxIev/o3oxRH97RY3DmsMxSNYeqSbThcXA2R0oBNFqfn51TnfnzLkAbp49tf+gTzvi7AgsM1N765ufCFNzbm9Wt5x8yoCS9+TEZOmgUAKJs/CZsyK4NhRmNxrNm4p3uC6bWnOkjP9jWUudoqoXKuryzbh0EXNpM5J21P9Deh2JxZikqnD71aJaFpjP30xHODsPHIrLtdLSLMZ3fSf8N5j1RJ4SYse3W+K/yG4ctVTbuQIZSoN7hjU1BKsHxXdm2c1cm7yjjPyHP54zl4Xt3vOt87A3srArRpurXn0Wrvs0Fdv1DXiTnGZsK9wzvDH1QxZckOqIxf0P6et2PynLqz/2Pvi5f3SPVLosgnXjG43m/4l09fh9XHpjGLJHziU+hVCuMd3l+xF8O7NsegTmmYdFl3PLbg7JwWAAACVDh9+GbDETw8qicm3zQAQzqmwhVQ0bV5PLqnx2N/bjm2HC2C06/g1/156Nw8Hi/cMACdmychoOro2yY5VB4VBAbp9OTEPm2TseKZ0circOGXPVn4Yt1hVPuUPwyF4hywyCImXdkbLZJisGjzESzaEQqdMojC0jiL9HNE82iMuKhDvfcvABxyBNF/dOuYbQWO2wIa6y0Q3OdWtMwCh39jnNX0Xd8HZ+3YMP0n3/FqN8LN0qVelb1Q4KY7c6qFV53B8APVC14CABT5GCyy0ETnvO2JsgY05DTbdbwEQ7o0x7U9M/DKoi2hsDAgqOlso+na1+FfOAnktVvq7ZrO26iWvToe1jFTIIBtI4T7OIM5MdKGnhlJcHj8+G7LUbDaqV8djPNEXVc6lXgCJ4wqt0ZHhdfZpdKrzVN01sJuktE3IwFX9MhA29RYrDuYD1kUoOp669xq/1y/zsUqn8pyK49OUutPz/F3JFo17J7D86OuFt9wKtq7eVVu+4sLN+GDu4fhhoEdsPxAHpbuyf1Dldo/hyCoali48TDG9m+L5glRaJ5wMpeHMR3fbT6KzJIacHDM+Gkn7CYZY/q2xgOXXwAA4EwHoYDT48eRwgqAcBRUOHE4vwKEcFhNMgZ1TMOQzs3QLS0BT3yxDiVO/++WgJRwjO7VEsO6NkdWcSVe+mYTXH4Fskhzo4z0tf15ZW6++AWQDx9okP4td3hhlmWTQHBEJCSBgTdRdNaNgHQrdAfGlroDu+xj+6+zSuJ6Z1C7VGUsg3Me7dfp9FOdRHk1NYg2GQcyzpvU/ZYA8PoVLNmWiUGdmmFgx6Z4e9nu2rpppMpAyTHZINV7tH29KIXYTRIowXGXolWAIzUhwowoqxFGWcKl3dJxtLQGVd6TLmLOYfIElcsfvjDmxz0tZum/Tr4LJpGAgjs5GOUIuXbdfhWvLtqGGl8Q/qAGhTFwwOxR2QgA8DMdleBrvWXufdc+Px8Lnx9b7zd997uPodvE92GVhW+2FdYM83M2dsnOLAzccAh3D++OV6/rj+OlDhwrqwEhwpkdlIce+glDOyIx0obcsmpsOFIIl19Fx9Ro9GyRhLEXtcNPu7Ox/XgZiqq9eGjeKvywLQv92jWBQAjapkRhWNd07M0rx/bjZYAgYsXeXAzO/gqMh4qfDe6QiklX98G4iztiT245Zvy887QAZ844ujaNwTOj+0KkBLOW7sLu3AoIlMAkCp/mzX1o25XPfgbSEIIytTx+1SW47bKJ+T3uvuKOSq/WrDoQuCqg0wGM8Y4a49E6+BAlyIf4Fd3NQWRCAJnytWkWHPJoQA2AYU/PweC0SPr0iqO9Ga97pgnAGZrG2NCvbQooBRIjLIi0GuEKKKCU5EZZjSWUENR3aHa9SHfEWQxItJpLBEIzCSHYnVeJZxesg9sXwKNX9sbsO4aF1JVOKWQQ4KTnV/t9sUcqQsoCvVMMyHuyQ5ZJpIspAQodHqzLLGK5Fc6g0xdw64yVE/B8k0i2R5gMCwVK3rbIwjNWg7jvhm7JtMbRMBvBALDj7Tuxp6Q6EG4gbxooyVR14I3F27FmXzY6pyfi2at6I8JoPPONVs6RkRiJa3q3htMXwIT3l+Kmt3/BPXNW4Mop3+P7LUeQnhCFS7o0AxCyA5/K8OOeHDz+6a/4dW82WidHQ9V0fLp6X0gFmAAKYyhyeFDi9CC70oX3V+7F4q2hqh4tkqJPV3ViHHFhJjw7ui+axkbgoxV7MXvlPoByGEXyc6SRzUoaNwXfv1T/gbOnctslrQBtBbbNvEd77cpemY75j7/WMsJ6RYxZHmw3iPeaRLKUEOrROaw6ZwYCrsmi9N3mMsW7f84jAID9xR7M2FwQpzPajddJOesMvVok4osHL8fofu1wILcMD89dieJqdyiSgmPfwffuccXY6l9OqV6M6uL2TbHt46V+gdLdIEBA1TH71/24e/YyHCupxpW9W+O92weja7M4cBYatDnnqZ6A2tZRq7D6zcv3wPzEZm6TyPxwo/SUURDvtUjSuCiL4XKjSIcl2QwD7Sb0T7aT4esev+z6cT2aPGoWsNYbUPuvKvP2XHG4Ens313tC8AkuTDCg+NOf91gk8UWJwpNb7sTzCzeguMqFMRe2xUOXdoNEyZlvIp4Q2AyFXnFdB/TTkwhPFnsJORPMkoCbB7TH9FsHIyUmDHNX7sW327JCSUE8VMgcnAOMgzOOJtF2tEkJTSur3f4TUfYcgCRSPHpZD4zo1gK7jxdhxs+h0B2J0iKbLL5a4ualxZ8+3mD9+UeMHtEDALB31t2+wrkP7HNs2DmrSYTlOpssXkMJyQIIKEixWSbbYu0nR0+nV4HXr7bROW/KEXqoL+3WDLMmDEWvVslYfzAfd81ejp/35kFjAAXhBlHYSS5/Fasn313v11Ev4/qbd1wC+ZrJkAW6NaCxIAM36Az4ZusxFFd7MPO2IRjUqRnmRdhw74fLse5oMXTOLarO+nq/XLby5mc+xycv3QD/N0+hANiF0A8AnFDrKR52N+KbdABHIKz/lEUDFI1dH1D1KxhHhCSwpMHt7XvGf77Rf25X8PcsmfoAOk2cBYNAvz9U5rpIB25Zf7SYvv7dJrx6w0W4+5Ju2J1XhkXbs39nWOQU5QgSUqJEdrkTq/fl4qo+rTHjlkH4csMhOHwKLmqTgst6ZqCoyonVBwsAkBN7SPcP74w7h3eFSAm+WLMfT3+5Do7a/ai2SZG4sX87+AIKKjxBRFoMGNapKfq0aYJjRRX4aedJxWAK4LKuzXDr4E4oqfFg0mdrkFnihCgQn0USpj3YJ3HzsmPVWNVQnXkG1C0VytxdHJJAj3PACACiQHY3izHmBTWgHMBtLy7EnGffh3XMsL6huuMUV/VogTfHD0RsmAW/7MjEvR/9WrtPGtL0FwQhP9wobLFIRjREJGe9TZY7JtggErpxd6lrb1DnPULPDsWmY6UYM30JZt46CMO7pOPOwR2xLasYPh1QGb8y9daRH60oLP/T8r19734PD6S1xVNHd4U5Vd/wcp9yp8bQTeewAIBA4KSc+71+YqAcDWZUALDn7bsQdcPrvjir9EyhU09XGC6avXIf0uLseOCynph562A4PD9i1YGC03buGWfwBVXojMEVDG0YO3xBPPzpKnAAw7s2xwvX9QcA6DrDsaJqPPrZWqw/Eqr4xxnHqO7N8dAVF6CsxoMp323Gh6sPwhMMefQ44+iUGo17L+kKi7FO/50jqKrYfKQAT85fi63Hy058tm/LJEy/ZTAMEsVTX2zEigOFEChgEOnXGfER732y18HqWwX4VOZvcmPs5WPQ/sphbQKM+9tE2PJcQYWtfuuO0z6XVaag/cTpsIisOQU3EwKfTZY+X7+9JPjJwyOwbRqwNLsMTcZfmVLm9V3JuI7kCBueH90XCRF2zF62E08vWBcSzax1lHEQiBRrn72kw+GNOTX4uAGur96MKkwQsfKDlWWmS7uuJwQ96pYXhBJkldRgV1YJhndJx5GS6pCICaXQGFo7A9qFDp86/7fHu/6FjxFnNdDvDpVm3Lxr82Wqjks1ji6cw0wIYTLFUUrIj0YRS9pFS9vX7i32fzy6G+rPMfrH9EmKwJL9JaUJUZbna4Lsw6DOW0z+fjsizCbcPKgTpt58Me7/eCXWHS48cSNVnWPaDzuwYm8elh8sPFEKNK/KgwkfLEOXtDi0To6CURSQX+XGzuxS5FS4To54hOCXPbkIvP8LdmWXYW9eJVR+stAbocDiXTmoeuM7NI0NR7hJhiugIqfMgZ3ZZahwh7Y0GAvtR70x7iLERVjxysL1mLtqPwACg0iWxZilF7KKKgJVX05q0D58cs5ctLxyeMJxp+8tjSO13Kt8YTXi68dfnnfoqFvniyaH7mJ6nIwmN0+DRLBOENmlOkMrgyCtTkmLxc0jQ17QIq8HYQahj8ZZa0Ioyl0BHC6sRIukKGw6UogKlx9EOOlAIgAXCTbe+tkWHYufapDrq9dVmnjNZJhFeoVX079iHCfCmQVK8Nm9l+Davm0wdtoSfLnxMIgQyng1isK8m7o3u/VgSY2+/s3bQl8Y+QySo8OaeYLamICOW1SGZiw0I1BFisNGQfgxzmL4uHOaPSe/MsA6xYjwqzotqlHEkjK3cmD+E/XuJj2V0c/MxlfvbEH8yFZXVwXYHFXnYanRVnwycSQuatcUO44VYdzMH3GopOYUjycHGH4XgcFDxStq1zsnK46fVqit9vucsdAeC/2DW8dr2+B16Yj0lGPVqrRGWfHR3Zfgovap+HzNAdz/ySo4fQpEgeRGmsTRLr+yzffVpAbtOwDA1ZMRaxUuq/bpcxnjkQCBJCDLINJPjIJwUAfxugKKP9wg+e0m2VvlVwK+YCCQFmGudvkVpeTTUOpJn/vfR+8Emc7cXf6hqpNbOEJryaeu6oWXr78QL3+9Ac98ufE0jT+BkIoIgzBUY3x3zRePNsjl1atwb7hMYRDIbgJyYjoXco8b0CIhAi5/MFRsmpxcgutMv3BtVmmrnErXieMQoxEKRzePhucVhuYgRDEK2BZpFG9PtBiGvTX2wqfT4mzHw1TOtmw8hv0lQcuPmY7bt5b4X9GtNlPrO99qkM6q46uXJqD9mHbolGhaFGWWX5QE6s2r8uCO93/B1iMF6NYiER/fOxztkqJOeDwJIX8Y0kRIbZkdQQi9aITf61WEiubVfubP9sMIao9Da48Vai8kDceRGG7FrAnDMLBjGr7ffBgPzVsNpy8IifKSKJN0b+msq7a1TbQ3vEHVnqtRoCujjeL1JoGuoGABVefp7qD+cpVf/brGF/xRY3x5TUBdme/wrnYF1PUEZJHOtGYiPWVvyuXBD7nuDM7R/8SeFec4WlQJzjlaJUVDOiUrmAAwULIlxWTMjDedfZmkM6VejSo1QsIVrU2FFokuOnFrOEer+DA0jQlDaZUbRdXuU7LBOTSGtEKH98bCTx7CkEkfAQCahhsQaRLWmUWy2CSSr+0GckOLaHlEy4hD86y0oGT8iJ58X0ENNpa4LBEtEkftLvN+4wjyaX5dv73G77/4SA7Q7e5XGqzTAGD/rIdwvEzTmltM75lEvCVSEswsdeK+j3/FugN56NkyBbMmDEGflkl/mvvU0BCERre0GBumjhuAwZ2b4esNh/DoZ+tQ5QlAFEiVUaIvXJJh/aXrI99jx4yGW0edxtePY1z7CF/JvIeWpUcarg03iTcaRPojAdE5hygJtEQgJB+AoDEeD86TDaKw44JES3bT8FBI0aVPf4DCuY+g2BW4SWNofvKiCbLLnHB4A0iPD0eYWT6x1UFAFIuBLtiVleftmlQ/arR/xBnuVp4ZJduW4ljcAG4QqTeg86s4YATnuPqClhjVqxW2ZhZi3tqDp2SkkzrtBHvSBZcuLvMGPe7dS+HYtRzLp7/lOVBU8WubuIhvi11+g1clLTrHtjwe4Hae2nuEVOHzd6706s95Vf0ZlfE2BCCSQI9JRN84KWbmkU2eDDj2b2ywjgOA6j3LkHThJVpqlGFbpVsL0zjvVljppvtyy3BBiwT0aJmCTmkx2JdbhuIaD+q9cNRfQAgH40BylB3Txg3E1b1b4ecdx/DgJ78it9IFUaA+u1F87qkrOszaU+JjG968q0HPh69V8a4ahrjuw9LDOw1NWfK6r2yvsgNrp98fiO4+8lCMxbDaF9Ra6uAtjSLdGmOS7rcY6Pc64/tEgRwNM5vmHigPFmbWRnZorYYisuvweIdfe17jiD+1LZ2xkKpXtB0r9uQgv9Jdq7tOcyLM0qv2CJtjw4xbG+xa671uQ1K4ETE2w36BhkKHKKFolxoLgCDcakJSpLU2A/YkGuetqwNq/yL3SZmouxYtg2Dwu37NKb+v2q/95PCr07cUVzatCaophyq8rzn8+Nmv6eMZh0WkpNgqS8/EmQ2XVHz++Lf7k5/jOfMnN1inncqWqffgeIXmSYm0vWCV6AJBoPruvCqMm/kjth4tRJfmiVjwwOUY2j61ds3zT5wVB9M5msWEYe49l2DUBRlYtPkwbnt/KfKqPBAFwRtplF9uF2X/4JMVWWz5Cw2TI3UqLx5bisr5i1ET0G4v86lzm91GM7aVhO53ulCKSj8vjTJLD0qE7vIq+tByf/B1Qli2LLCp3i8eva9wzr27HJ8/eOJ4+dUeOD3+vqqut6l114BzDgqC5jHhMMoiImxm9GuVjLpOFwVsnHBBbEG7uLOs5HKW1OtIBQCDh1+H1ZtzgpZIU6Km84GEcMTZTejUNBZtU2MxoE0qNF3HwYJKaOyE2KQgEWZrFyUujeoyxFexcznS+16CVW/nMmsbW6ugzq5gHDF+jber8eu3+FR2BQO3CoRUmSRxboRJfLh9XPTiIIPz3kt74NCaJQ3aab/l+bcW47Mfv/G3TojZ4FX0OI3zDiVOH9mWWYTkSCt6tkxBv9YpUIIq9uSWQz/Fc9cQcMbQs3kC3r71YlzYLhWzl+/G01+sQ6nDB0mkHrssvdYlKW5abo03mDnn3n+kjy68bAzW+s0wW+XYgMZu8WmsdYTJsrsmvk/FV/cMRn65A1v35tdERpnzFYY+GkNnVSOdoszy/gfnri6pyj2MEf26AQBa3zkNCTYxqszHXlYZaYnaaW5cmAUPj+iKl8deiGbxUTiQV4r56w/iWJkTIqU1kUbpxV8zHZmH3q/PKsW/p96N6sC67xE98DKYDYLHr/KhjCF8V0459ueWIz7MhG4tEjCoY1M0ibQhv9KJCpc/5LDgaBLQyf78quCBqfdPxOynbkar8aMRYRRyHAG9j8pJms7RTOdIoAQ+o0hX2A3i05d1SH0vym4uXvrKjaxyx48N2ll/xrK5LwHHNsHUcagvOcywxRPUw1SGjqUOL916rBgJYRZ0z0jChW1SIAkUB/Mr4AtqDWBYHCKhGNohFTNuHYz2TWPx7s878OI3m1Hq9kMSqNsmi68NTo6Z4Q0GA7tmNeyU71TWLvkSrQeNQpRJynEpekxQ49cENK1jbLghc0e2q7jYGeS3DUjB0tf3Z4V1ivQojA/UGDI4eGXg68OrBl7ZDht+/gaZe/bjuZ8OQgO9zKvwhziHIIsE/Vol4c1xAzH+4k6QRIqvNx7Go5+tqd3ro5AFrE61C2/YDESt2Lm8Qa+13o0KAO685x4sf/GWUnOHwckKQ28Q4Hi5Az9sz0KFw4XOzeIxoEMaRnRuBn8wiJ055eCAwDlMFzUNX7yzxKEWb/8F48fejO/3lgXCrAYe1NlIDggSJSVRJumu9GjrK8dKCg52T26uL3j66gbtpDPFsWsp9NYXe2PthvW+oGJmoB1rvIq0bO9xmEQBF7RMwsAOTZERF4Z9eWWocAVqIyzOv23OOCyyiLsGd8Bbtw+B3STj6S/WYcrirfAqGiRCqu1G4ZE+aTHvZ7sDwc0z7zj/Rs+Syp1L4W01WImzGdd7FTVW1fllfk27vCagqmaTdOT7Qh5o2zMcaRHGgxU+RQDhLewG4XWhTVT+7VcPxTcfvoVtUiu0izGZj1erL6o6a2OQBTx7ZW9MGTcQHZrGY29OKR74aAWm/7wdhdXekEcUYAaBTCso92/5+IZOWPDF/PO/mL+gQYxq8w9fwjR6MkwiVYI6H8UBIyEEAU3HjuxybDxciAiLERe0SoEkivhm81GooXVWolvh+Zl7y/dMfe1JPHfXKKT1vRzRZkOeO6g21xhrRwh0WRKWZBW69i588go8cf2ABu2gs8W7dznMbQYHwwVlPRckr8Z414DKTRuPFqGgwoluzePRo2UK+rZOhsvrx9GiqvOcDobUlZrHh2PKjRfhvpE9kF1cjYkfrcBXm49C44AsCtlWmT4cydUFhQ6fmvnRA/9a/yj7lsHSZlDQKgkbOEW1wni/oM5HuAJKPxsUg90olmz/bp0joW3T/WaBbmkRKWxOtIn6m+OH480Pf8Hcn/fBJcljvRq/h3Eux9rNmHzjRYgLs2Lmj1vw+OfrsCmzGDo/udcnEXI80iy9aDHJNR8+3rABwkADGRUAtO53OWItpipnUO2g6rwNENqr4QAKK1zwBTVc26cVDhdU4KtNR2srUxBZ5zyqWVr4D5tyKn3O3Uvh2L0Ujqa9AzajVBbU2EidkygA5OKM6B/WHyzW8jf9O1O+v8K5dxkGXXWzOrRV/Lb9RdVlOkfngMrC9uZWYH9uGTISwtGleSIGdWwKi0yRU+5EjefMCgecSl2Q7YguaXjrlkEY2CENy3Ydx8SPVmBDZjE4IdwgkF1hBuGhys8fXfzCo5PYVy82/EP1d7j2LUOXIdcEWkZat1X6lN0qY7EaQ09FZ1c4A/oAQ9OkVqrGfJc1C1sFatB+mRxy9R+K7oLIOHtUlU+bonGeHkrQpLi+bxsYJQEPfvwrDhdVgZ4eQcGsBnHmR1devCTL4eK5Gxc1+PU1RNVWAMBt/dOxLafaY5Okdykhp+dlEIIIqwxRoKjyBKGyuuhpDoXxnjV+dXje7iK89M5CAEDXlDD0SbZuNkvkY7NEv4wzm14LM9KgTBvs9M+bhc/ehN0lTm3mzT0+ibXJN8gi2csJsHx/Lq6bthifrNwNSRQx6Zp++Oy+kRjULgUypbWKvn+nJxGKRE+OtOKVMX3x+QOXIyMxEq9/sxHj3vkZBwurQCllBpGvSLDS68q3ly6/8ZEZmHhT/WdHnyur37gNN/VM1sqzy5Y3DxeujTEKVxkF8WuNw64zfp1FFJLmTNnJ1NotvhFPfoDCuZ/CG9Qu1xgPhbMTwK9oqPEpMEgSTCYDQChOTVIUCPKMMv1i/KKVbPWUcf/ItTXoxkmL29+FTKkly+FaqOj8krrfc13HHYM74f27hmP6D1vx0CdrTgslMQhkX7zZMCagscOlnz0EAIi66VU0izBZgyphuS6vz/nZk//M7v95MuH1zzB70jQk3ji2tUvTHw2ouEHjTLIZJIzo0gyTruyFjmnxqHZ78f2Wo3h32W7szimvq1h42rE4B8AYYsPMuK5Pa9w2uCNaJkVhw8F8vPLdZmw4UoSgziBR6jHLdHqEWZ6V6+MlVxor8N1HDbsZfj6sWnsAA/v3w5DHp0uHKxzRMiF2E9Q8xnng0MehUKKYG6bBJJK25b7g9wpDizrDkUUB3z18OQZ1TMMlr3yNX/flnghmJgCMAmY+clHCgz8ecbFdb9W//v4f0aBPJecBSGNmwCZLV7mD2oc64xG1f8Br1/bG49deiJe/WotnFm46Ld2egsMiCVOe6RT5xJfH3Wzn7AfP9RT+MzS/dSosBjG8wBm806PyiRrjiZxxtEuOwgMju+KqXq0RbjXhcH4FZq/Yje+2ZqKg2o26yQTnHFE2Iwa3b4KbB7THoE7NUe70YO7KvZi1bDcKq32gAuUyJYfsBjqlc5J1YY1PCWybefbu49Rx0yBQYghy0rTSp0SaROSlR1uqMiu9F4OxQK84cU21orHtsxo28LaOjrdPxQ0ZBvrSTt+rHo09fupWn0AoPps4HKP7tsG1by7CN5uPnjAqiZLMGLNwpV/VD1bP/+dywxp0/kSIEX2S7BiRbvvBSPFLnQVLooDUuHAAHBXuAKAznCpzw0EQ0HHTe5nui3eWe/HYtK/+sQ75LX3veQ+XPvmpPOzhD8WhD394zsc5/tHDiDWbHU/1TZ8cbzaMMwlkIxUIO1hYiYkfrcTNM3/AjmNFaJEUiTfHD8Kix6/CjRe2Q4zdhFi7CaMvSMfXD47EvPsvQ/92TbFyTxaumvw9nlu4CUUOP0SBqCaRfJ1ok0aVffrwp0KAnZNBAUBQU+FXgilVvsBXKuOrGMij+U7lgqDO56scs3KDYrOyYIMtx0/jjZnfYm81w7uH1H5+xm881aA459A0DQ6vAkIorKfE8xEAIsHCoq+PH7yoWcMUVPgzGrxncjf9iDHvLtdtsuAO6vxyBhiMEsXtgzshLS4CFoMIRWPIrXAiqKgg"
               +
              "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",
          fileName=
              "modelica://TJUthermo/../../倪先生的毕业设计/天大图像/u=1379584405,4084835147&fm=21&gp=0.jpg"),
        Rectangle(
          extent={{-76,42},{8,38}},
          lineColor={28,108,200},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={255,255,255}),
        Rectangle(
          extent={{-76,24},{8,20}},
          lineColor={28,108,200},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={255,255,255}),
        Rectangle(
          extent={{-74,8},{10,4}},
          lineColor={28,108,200},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={255,255,255}),
        Rectangle(
          extent={{-74,-10},{26,-14}},
          lineColor={28,108,200},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={255,255,255}),
        Text(
          extent={{-70,-22},{70,-70}},
          lineColor={215,215,215},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={255,255,255},
          textString="TJUthermo")}));

end TJUthermo;
