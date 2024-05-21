classdef TransmisorFMT < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    % Atributos del transmisor FMT
    properties
        potencia;
        pulso_conformador; 
        parametros;
        
    end
    
    % Señales de salida de cada etapa del transmisor
    properties
        sk;
        impulsos;
        xbb;
        xi;
    end

    methods
        function obj = TransmisorFMT(parametros,potencia,pulso)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.parametros = parametros;
            obj.potencia = potencia;            
            obj.pulso_conformador = pulso;
            fi = obj.calcular_portadoras_bb(obj.parametros.Bw,obj.parametros.numTonos);
            obj.parametros.setFi(fi);
        end

        function x = transmitirBits(obj,bits)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.asignarPotencia();
            obj.sk = obj.mapeo(bits,keys(obj.parametros.constelacion),obj.parametros.M);
            simbolos = keys(obj.parametros.constelacion);
            obj.sk = obj.conversor_serie_paralelo(obj.sk,obj.parametros.N_relativo,simbolos(1));
            obj.sk = obj.parametros.distribuidor.atenuacion(obj.parametros.distribuidor.atenuacion>0).*obj.sk.*sqrt(obj.parametros.N_relativo);
            obj.impulsos = upsample(obj.sk,obj.parametros.upsampleFactor*obj.parametros.numTonos);
            obj.xbb = obj.filtro_conformador(obj.impulsos,obj.pulso_conformador);
            obj.xi = obj.modulador(obj.xbb,obj.parametros.fi,obj.parametros.frecMuestreo);
            x = sum(obj.xi,2);
        end

        function x = transmitirBitsClear(obj,bits)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.asignarPotencia();
            sk = obj.mapeo(bits,keys(obj.parametros.constelacion),obj.parametros.M);
            clear bits
            simbolos = keys(obj.parametros.constelacion);
            sk = obj.conversor_serie_paralelo(sk,obj.parametros.N_relativo,simbolos(1));
            sk = obj.parametros.distribuidor.atenuacion(obj.parametros.distribuidor.atenuacion>0).*sk.*sqrt(obj.parametros.numTonos);
            impulsos = upsample(sk,obj.parametros.upsampleFactor*obj.parametros.numTonos);
            clear sk
            xbb = obj.filtro_conformador(impulsos,obj.pulso_conformador);
            clear impulsos
            xi = obj.modulador(xbb,obj.parametros.fi,obj.parametros.frecMuestreo);
            clear xbb
            x = sum(xi,2);
        end

        function asignarPotencia(obj)
            obj.parametros.setConstelacion(keys(obj.parametros.constelacion).*sqrt(obj.potencia.*obj.parametros.T));
        end
        
        function simbolos = mapeo(~,bits,constelacion,M)
            n = log2(M); %numero de bits agrupados por simbolo
            if mod(size(bits,1),n) ~= 0
                error("Longitud = %d, La longitud de los bits de entrada debe ser múltiplo de %d",length(bits),n)
            end
            
            indices = bit2int(bits,n) + 1; %Obtiene los valores decimales de los grupos de n bits
                                           % y los escala a valores entre 1 y 2^n
            simbolos = constelacion(indices); % Extrae el valor de la constelacion ubicado
                                              % en el indice especificado, debido
                                              % al orden ascendente de la
                                              % constelacion
        end 
        
        function out = conversor_serie_paralelo(~,sk,N,simbolo)
            num_sk_extra = N-mod(length(sk),N);
            num_sk_extra(num_sk_extra == N) = 0;
            sk_extra = ones(num_sk_extra,1).*simbolo;  % Crea simbolos extra para tener una cantidad 
                                                      % de simbolos multiplo de el numero de tonos
            out = reshape([sk;sk_extra],N,[]).';
        end
        
        function out = filtro_conformador(~,in,p)
            M = size(in,1);
            N = size(p,1);
            L = M + N - 1;
            out = ifft(fft(in,L,1).*fft(p,L,1),[],1);
        end

        function fi = calcular_portadoras_bb(~,Bw,numTonos)
            fi = (Bw/numTonos)*((-numTonos/2:numTonos/2-1)+(1/2));
        end

        function x = modulador(~,in,fi,Fs)
            L = size(in,1); % Longitud de la señal
            arg = ((1:L).'*1i*2*pi/Fs).*fi;
            x = exp(arg).*in;
        end
    end
end