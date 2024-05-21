try
    delete(loadbar)
catch
    pause(0.000000001)
end
%% Inicializando parámetros
a = 1;
constelacion = [-a;a];
M = length(constelacion);
T = 1;
upsampleFactor = 2;
frecPortadora = 1.5;
rolloff = 1;
numVecinos = 1; % Numero promedio de vecinos de la constelación
ro = 4; % Constante que relaciona distancia minima con Eb

N = 32;
Bw = (1+rolloff)/T;

alpha = [1,0.5];
% alpha = [1,0.8];

% Tau = [0,T];
% Tau = [0,3*T/2];
% Tau = [0,2*T];
Tau = [0,5*T/2];

parametros = Parametros(constelacion,T,upsampleFactor,frecPortadora,N,Bw,alpha,Tau,numVecinos,ro);

%% Parametros de pruebas
potencia = 1; % potencia del transmisor
EbNo_dB = 0:11;
EbNo = 10.^(EbNo_dB/10);
No = (potencia*T/log2(parametros.M))./EbNo;

%% Canal
ch = Canal(parametros,0); 

%% Transmisor


periodos = 64;
pulso = rcosdesign(rolloff,periodos,upsampleFactor*N)'; % genera el pulso de raiz de coseno alzado de energía 1;

transiente = periodos*upsampleFactor*N;
parametros.transiente = [transiente,transiente];

tx = TransmisorFMT(parametros,potencia,pulso);
parametros.setHi(ch.h)



%% Receptor

rx = ReceptorFMT(parametros,pulso);

%% Transmision
BER_teorica = qfunc(sqrt(2*EbNo));
% semilogy(EbNo_dB,BER_teorica)

k = (EbNo_dB+1).^3;

numbits = 50*round((1./BER_teorica))./k;
n = 2.^nextpow2(numbits/(N*log2(parametros.M))+2*periodos);
numbits = (n-2*periodos)*log2(parametros.M)*N;
BER = ones(size(BER_teorica));
porcentaje = 0;
bits_enviados = 0;
bits_totales = sum(numbits.*k);

loadbar = waitbar(0,'Inicializando');
global stop 
stop = false;
loadbar.CloseRequestFcn = @(src,event)closed(src,event,data);

letters = ['A':'Z', '0':'9', 'a':'z'];
porcentaje_tot = 0;
for Tau = 2:5
    parametros.Tau = [0,Tau*T/2];
    ch = Canal(parametros,0); 
    parametros.setHi(ch.h)
    name = strcat("test_05_",num2str(Tau),"T2_",letters(randi(numel(letters),[1,6])),".csv");
    for metodo = 1:6        
        for i = 1:length(EbNo_dB)
            BER_aux = zeros(k(i),1);
            ch.No = No(i);
            parametros.distribuidor.seleccionar_distribucion(metodo,parametros,ch,tx)
            
            for j = 1:k(i) 
                bits = randi([0 1],numbits(i),1);
                
                x = tx.transmitirBitsClear(bits);
    
                % Canal
                y = ch.aplicar_canal(x);
    
    
                bits_recibidos = rx.decodificarBitsClear(y);
                bits_recibidos = bits_recibidos(1:length(bits));
                BER_aux(j) = sum(bits_recibidos ~= bits)/length(bits);
                porcentaje_tot = porcentaje_tot + numbits(i)/(6*4*bits_totales);
                if stop
                    return
                end
                waitbar(porcentaje_tot,loadbar,strcat("Transmitiendo bits, Escenario ",num2str((Tau-2)*6+metodo)," de 24"));
                pause(0.000001)
            end
             
            BER(i) = mean(BER_aux);
        end
        
        writematrix(BER,name,'WriteMode','append')
    end
end
% hold on
% semilogy(EbNo_dB,BER)

function closed(src,event,data)
    answer = questdlg(sprintf("¿Está seguro de cerrar el programa?\n ¡NO se guardara el progreso!"), ...
	'Advertencia', ...
	'CONTINUAR EJECUCIÓN','SALIR','CONTINUAR EJECUCIÓN');
% Handle response
    switch answer
        case 'SALIR'
            global stop
            delete(src)
            stop = true;
    end
end