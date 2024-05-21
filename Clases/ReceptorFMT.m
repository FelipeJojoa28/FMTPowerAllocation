classdef ReceptorFMT < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        pulso;
        parametros;
    end

    properties
        xi;
        xbb;
        uk;
        Uk;
        sk;
    end

    methods
        function obj = ReceptorFMT(parametros,pulso)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.parametros = parametros;
            obj.pulso = pulso;

        end

        function bits = decodificarBits(obj,senial)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.xi = obj.demodulador(senial,obj.parametros.fi,obj.parametros.frecMuestreo);
            obj.xbb = obj.filtro_acoplado(obj.xi,obj.pulso);
            obj.xbb = obj.xbb(obj.parametros.transiente(1)+1:end-obj.parametros.transiente(2),:);
            obj.uk = downsample(obj.xbb,obj.parametros.upsampleFactor*obj.parametros.numTonos);
            obj.Uk = obj.ecualizar(obj.uk,obj.parametros.Hi(obj.parametros.distribuidor.atenuacion>0));
            obj.Uk = obj.Uk./(obj.parametros.distribuidor.atenuacion(obj.parametros.distribuidor.atenuacion>0)*sqrt(obj.parametros.numTonos));
            obj.sk = obj.decision(obj.Uk,keys(obj.parametros.constelacion));
            obj.sk = reshape(obj.sk.',numel(obj.sk),1); % Conversor serie-paralel
            bits = obj.demapeo(obj.sk,obj.parametros.constelacion);
        end

        function bits = decodificarBitsClear(obj,senial)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            xi = obj.demodulador(senial,obj.parametros.fi,obj.parametros.frecMuestreo);
            clear senial
            xbb = obj.filtro_acoplado(xi,obj.pulso);
            clear xi
            xbb = xbb(obj.parametros.transiente(1)+1:end-obj.parametros.transiente(2),:);
            uk = downsample(xbb,obj.parametros.upsampleFactor*obj.parametros.numTonos);
            clear xbb
            Uk = obj.ecualizar(uk,obj.parametros.Hi(obj.parametros.distribuidor.atenuacion>0));
            clear uk
            Uk = Uk./(obj.parametros.distribuidor.atenuacion(obj.parametros.distribuidor.atenuacion>0)*sqrt(obj.parametros.numTonos));
            sk = obj.decision(Uk,keys(obj.parametros.constelacion));
            clear Uk
            sk = reshape(sk.',numel(sk),1); % Conversor paralelo-serie
            bits = obj.demapeo(sk,obj.parametros.constelacion);
        end

        function out = filtro_acoplado(obj,in,p)
            M = size(in,1);
            N = size(p,1);
            L = M + N - 1;
            out = ifft(fft(in,L,1).*fft(p,L,1),[],1);
        end

        function x = demodulador(~,in,fi,Fs)
            L = size(in,1); % Longitud de la señal
            arg = ((1:L).'*1i*2*pi/Fs).*fi;
            x = exp(-arg).*in;
        end

        function Sk = decision(~,Uk,constelacion)

            constelacion = repmat(constelacion,1,size(Uk,1)); %realiza replicas de la constelacion segun el numero de simbolos

            Sk = zeros(size(Uk)); % crea una matriz para precargar memoria
            for k=1:size(Uk,2) % Recorre los diferentes canales
                a= abs(constelacion-Uk(:,k).'); % Analiza la distancia entre todos los simbolos de la constelacion y el valor Uk
                % obtiene el valor de la magnitud
                [~,i] = min(a,[],1);    %Extrae el valor minimo en todas las magnitudes, es decir obtiene la magnitud mas pequeña entre Uk y un Sk

                Sk(:,k) = constelacion(i); % Compara el valor minimo obtenido y lo extrae de la constelacion, guardandolo en Sk
            end

        end

        function bits = demapeo(~,Sk,constelacion)
            n = log2(length(keys(constelacion))); % Numero de bits agrupados
            enteros = constelacion(Sk);
            bits = int2bit(enteros,n); % Convierte los valores decimales a valores binarios
        end

        function uk = ecualizar(~,uk,Hi)
            uk = uk./Hi;
        end
    end
end