clear
clc
conv_poly=[561,753];
conv_constraints_length=9;
trellis=poly2trellis(conv_constraints_length,conv_poly);
tail=zeros(1,conv_constraints_length-1);
tb_len=50;
N_chip=38400;
SF=64;
sf_code=zeros(1,SF);
sf_code(2:2:SF)=1;
sf_code_anti=1-2*sf_code;
Num_coded_bits=N_chip/SF;
N_bits=Num_coded_bits/2-conv_constraints_length+1;
Nuser=3;
SNR_db=-21:-13;
SNR=10.^(SNR_db/10);
N_block=10000;
ts=1/(3.84e6);
fd=100;
tau=[0 ts 2*ts];
path_position=[0,1,2];
pd=[0 0 0];
for nu=1:Nuser
channel(nu)=rayleighchan(ts,fd,tau,pd);
channel(nu).ResetBeforeFiltering=1;
channle(nu).StorePathGains=1;
end
source=zeros(Nuser,N_bits);
rake_fingers=zeros(Nuser,3,Num_coded_bits);
combined=zeros(Nuser,Num_coded_bits);
combined_real=zeros(Nuser,Num_coded_bits);
interleaved=zeros(Nuser,Num_coded_bits);

gold(1,:)=Gold(1234,N_chip);
gold(2,:)=Gold(76543,N_chip);
gold(3,:)=Gold(37657,N_chip);
for loop_snr=1:length(SNR_db)
err=0;
err_blk=0;
sigma=sqrt(1/SNR(loop_snr)/2);
chips=zeros(Nuser,N_chip);
pass_channel=zeros(Nuser,N_chip);
temp_chip=zeros(1,N_chip);
for loop_block=1:N_block
for nu=1:Nuser
source(nu,:)=randsrc(1,N_bits,[0 1]);
code=convenc([source(nu,:),tail],trellis);
interleaved(nu,:)=matintrlv(code,20,length(code)/20);
for i=1:length(interleaved(nu,:))
for ii=1:SF
temp_chip((i-1)*SF+ii)=xor(interleaved(nu,i),sf_code(ii));
end
end
chips(nu,:)=(1-2*temp_chip).*gold(nu,:);
chips1=chips(nu,:);
pass_channel(nu,:)=(filter(channel(nu),chips1)).';
end
received=sum(pass_channel,1)+(randn(1,N_chip)+1i*randn(1,N_chip))*sigma;
for nu=1:Nuser
for nsym=1:length(interleaved)
combined(nu,nsym)=0;
index_offset=(nsym-1)*SF+1;
temp_code_scrambling=gold(nu,index_offset:(index_offset+SF-1));
temp_code=conj(temp_code_scrambling).*sf_code_anti;
for n_finger=1:3
tmp_index=(index_offset:(index_offset+SF-1))+path_position(n_finger);
if(tmp_index(length(tmp_index))>length(received))
break;
end
temp_r=received(tmp_index);
temp_code=temp_code(1:length(temp_r));
rake_fingers(nu,nsym,n_finger)=sum(temp_r.*temp_code);
rake_fingers(nu,nsym,n_finger)=rake_fingers(nu,nsym,n_finger)*conj(channel(nu).PathGains(2,n_finger));
combined(nu,nsym)=combined(nu,nsym)+rake_fingers(nu,nsym,n_finger);
end
combined_real(nu,nsym)=real(combined(nu,nsym));
end
deinterleaved=matdeintrlv(combined_real(nu,:),20,length(code)/20);
decision=vitdec(deinterleaved,trellis,tb_len,'term','unquant');
decision=decision(1:size(source,2));
err_b=sum(decision~=source(nu,:));
err=err+err_b;
err_blk=err_blk+(err_b>0);
end
if(err_blk>=100)
break;
end
end
ber(loop_snr)=err/(N_bits*loop_block*Nuser);
bler(loop_snr)=err_blk/(loop_block*Nuser);
end
semilogy(SNR_db,ber,'-^',SNR_db,bler,'-o');
grid on;
xlabel('SNR(db)');
ylabel('BER &BLER');
legend('BER','BLER');
