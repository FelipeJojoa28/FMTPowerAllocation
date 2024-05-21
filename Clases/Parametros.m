classdef Parametros < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        constelacion;
        M;
        upsampleFactor;
        frecMuestreo;
        frecPortadora;
        T;   
        numTonos;
        Bw;
        transiente;
        alpha;
        Tau;
        fi;
        Hi;
        distribuidor;
        numVecinos;
        ro;
        N_relativo;
    end

    methods
        function obj = Parametros(constelacion,T,upsampleFactor,frecPortadora,N,Bw,alpha,Tau,numVecinos,ro)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            obj.constelacion = obj.normalizarConstelacion(constelacion);
            obj.upsampleFactor = upsampleFactor;
            obj.frecPortadora = frecPortadora;
            obj.T = T;
            obj.frecMuestreo = obj.upsampleFactor/obj.T;
            obj.M = length(keys(obj.constelacion));
            obj.numTonos = N;
            obj.Bw = Bw;
            obj.alpha = alpha;
            obj.Tau = Tau;
            obj.distribuidor = DistribuidorPotencia();
            obj.numVecinos = numVecinos;
            obj.ro = ro;
        end

        function out = normalizarConstelacion(~,in)
            E = sum((abs(in)).^2)/size(in,1);
            out = in.*sqrt(1./E);       
            out = dictionary(out,(0:length(out)-1).');
        end
        
        function setConstelacion(obj,constelacion)
            obj.constelacion = dictionary(constelacion,(0:length(constelacion)-1).');
        end
        
        function setFi(obj,fi)
            obj.fi = fi;
        end
        
        function setHi(obj,h)
            obj.Hi = obj.extraerGanancias(h,obj.fi,obj.frecMuestreo);
        end

        function Hi = extraerGanancias(~,h,fi,Fs)
            L = length(h)*5000;
            H = fftshift(fft(h,L));
            f = (Fs*(-L/2:L/2-1)/L)';
            [~, idx] = min(abs(f-fi),[],1);
            Hi = H(idx).';
        end

    end
end