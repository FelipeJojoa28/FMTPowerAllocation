classdef Canal < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
       No;
       h;
       parametros;
       awgn = NaN;
       z;
    end

    methods
        function obj = Canal(parametros,No,h)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.parametros = parametros;
            obj.No = No;
            if nargin > 2
                obj.h = h;
            elseif nargin == 2
                obj.h = obj.respuestaImpulsoCanalDispersivo(obj.parametros.alpha,obj.parametros.Tau,obj.parametros.frecPortadora,obj.parametros.frecMuestreo);
            end
        end

        function output = aplicar_canal(obj,input)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            output = obj.aplicar_respuesta_impulso(input,obj.h);
            if isnan(obj.awgn)
                output = obj.aplicar_awgn(output,obj.No);
            else
                output = output+obj.awgn;
            end
        end

        function output = aplicar_respuesta_impulso(~,input,h)
            M = size(input,1);
            N = size(h,1);
            L = M + N - 1;
            output = ifft(fft(input,L,1).*fft(h,L,1),[],1);
            output = output(1:M); 
        end

        function y = aplicar_awgn(~,x,No)
            sigma = sqrt(No/2); % Desviacion estandar
            obj.z = sigma.*randn(size(x),"like",1i)*sqrt(2);  % ruido AWGN
            y = x+obj.z; %Adiciona el ruido en la señal de entrada
        end

        function h = respuestaImpulsoCanalDispersivo(~,a,Tau,fc,fs)
            afc = a.*exp(-1i*2*pi*fc*Tau);
            idx = (Tau.*fs)+1;
            h = zeros(max(idx),1);
            h(idx) = afc;
        end
    end
end