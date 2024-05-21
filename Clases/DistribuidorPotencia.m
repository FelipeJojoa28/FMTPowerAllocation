classdef DistribuidorPotencia < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        atenuacion;
    end

    methods
        function obj = DistribuidorPotencia()
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.atenuacion = 0;
        end
        
        function distribucionUniforme(obj,numTonos)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            a = ones(1,numTonos);
            obj.atenuacion = a.*sqrt(1/numTonos);
        end

        function distribucionProporcionalDirecta(obj, Hi)
            a = abs(Hi)/sum(abs(Hi));
            obj.atenuacion = sqrt(a);      
        end

        function distribucionProporcionalInversa(obj, Hi)
            a = (1./abs(Hi))/sum(1./(abs(Hi)));
            obj.atenuacion = sqrt(a);      
        end

        function distribucionProporcionalDirectaCuadratica(obj, Hi)
            a = abs(Hi).^2/sum(abs(Hi).^2);
            obj.atenuacion = sqrt(a);      
        end

        function distribucionProporcionalInversaCuadratica(obj, Hi)
            a = (1./abs(Hi).^2)/sum(1./(abs(Hi).^2));
            obj.atenuacion = sqrt(a);      
        end
        
        % function distribucionWaterFilling(obj, Hi, No, potencia)
        %     % Metodo de bisección
        %     gn = abs(Hi).^2;
        %     lower = min(No./gn);
        %     upper = (potencia + sum(No./gn))/length(Hi);
        %     min_error = 10^-7;
        %     pn = zeros(size(Hi));
        %     while abs(potencia - sum(pn)) > min_error
        %         lambda = (upper + lower)/2;
        %         pn = lambda - (No./gn);
        %         pn(pn<0) = 0;
        %         if potencia < sum(pn)
        %             upper = lambda;
        %         else
        %             lower = lambda;
        %         end
        %     end
        %     a = pn/potencia;
        %     obj.atenuacion = sqrt(a);
        % end

        function distribucionMinimaBER(obj, Hi, No, numVecinos, ro, potencia, M, numTonos, T)
            alpha = numVecinos/log2(M);
            beta = (abs(Hi).^2)*ro*potencia*numTonos*T/(2*No*log2(M));
            [betaN, idx] = min(beta);
            alphabeta = alpha*beta;
        % Metodo de bisección
           
            lower = 0;
            upper = 1;
            min_error = 10^-9;
            g = 20; % Funcion objetivo
            i = 0;
            exponent = 10000;
            

            while i < 2
                pN = (upper + lower)/2;
                g = sum(lambertw((((alphabeta/alphabeta(idx)).^2 * betaN*pN).^(1/exponent)*exp(betaN*pN/exponent)).^exponent)./beta) - 1;
                
                if g > 0
                    if g == Inf
                        i = 0;
                    end
                    upper = pN;
                else
                    lower = pN;
                end
                i = i+1;
            end
            % metodo Newton Raphson
            while abs(g) > min_error
                W = lambertw((((alphabeta/alphabeta(idx)).^2 * betaN*pN).^(1/exponent)*exp(betaN*pN/exponent)).^exponent);
                g = sum(W./beta) - 1;
                g_1 = (1+betaN*pN)*sum(W./(beta*pN.*(1+W)));
                pN = pN - g/g_1;
            end
            a = lambertw((alphabeta/alphabeta(idx)).^2 * betaN*pN*exp(betaN*pN))./beta;
            obj.atenuacion = sqrt(a);
        end

        function seleccionar_distribucion(obj,i,parametros,ch,tx)
            switch i
                case 1
                    obj.distribucionUniforme(parametros.numTonos);
                case 2
                    obj.distribucionProporcionalDirecta(parametros.Hi)
                case 3
                    obj.distribucionProporcionalInversa(parametros.Hi)
                case 4
                    obj.distribucionProporcionalDirectaCuadratica(parametros.Hi)
                case 5
                    obj.distribucionProporcionalInversaCuadratica(parametros.Hi)
                case 6
                    obj.distribucionMinimaBER(parametros.Hi, ch.No, parametros.numVecinos, parametros.ro, tx.potencia, parametros.M, parametros.numTonos, parametros.T)
                otherwise
                    disp("Error, no existe ese método de distribución")
            end
            parametros.N_relativo = length(obj.atenuacion(obj.atenuacion > 0));
            fi = tx.calcular_portadoras_bb(parametros.Bw,parametros.numTonos);
            parametros.setFi(fi(obj.atenuacion > 0));
        end

    end
end