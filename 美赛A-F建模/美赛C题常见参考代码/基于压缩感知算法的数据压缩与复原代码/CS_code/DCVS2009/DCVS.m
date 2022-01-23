function DCVS(SeqName, NF, GOPSize, MKP, MCSP, FrameWidth, FrameHeight, Iteration)
Ave_PSNR = 0.0;
Ave_K_PSNR = 0.0;
Ave_CS_PSNR = 0.0;
Ave_num_nz_theta = 0.0;
Ave_num_iteration = 0.0;
KFN = 0;
CSFN = 0;
MK = MKP * FrameWidth * FrameHeight;
MCS = MCSP * FrameWidth * FrameHeight;

tic; % start time

for fn = 1: NF
    if mod(fn - 1, GOPSize) == 0 %key frame
        if fn - 1 <= 9
            imgfile = sprintf('%s_00000%d.bmp', SeqName, fn - 1);
            Recimgfile = sprintf('Rec%s00%d.bmp', SeqName, fn - 1);
        elseif (fn - 1 >= 10) && (fn - 1 <= 99)
            imgfile = sprintf('%s_0000%d.bmp', SeqName, fn - 1); 
            Recimgfile = sprintf('Rec%s0%d.bmp', SeqName, fn - 1); 
        else
            imgfile = sprintf('%s_000%d.bmp', SeqName, fn - 1);
            Recimgfile = sprintf('Rec%s%d.bmp', SeqName, fn - 1);
        end            
        [Rec_K, K_PSNR, num_nz_theta, Num_Iteration] = fast_cs2d(imgfile, MK, 'BWHT', 0, 32, Iteration);
        Ave_PSNR = Ave_PSNR + K_PSNR;
        Ave_K_PSNR = Ave_K_PSNR + K_PSNR;
        Ave_num_nz_theta = Ave_num_nz_theta + num_nz_theta;
        Ave_num_iteration = Ave_num_iteration + Num_Iteration;
        imwrite(Rec_K, Recimgfile, 'bmp');
        KFN = KFN + 1;
    end
end

for fn = 1: NF
    if mod(fn - 1, GOPSize) ~= 0 % CS frame
        if fn - 1 <= 9
            imgfile = sprintf('%s_00000%d.bmp', SeqName, fn - 1);
            Recimgfile = sprintf('Rec%s00%d.bmp', SeqName, fn - 1);
            SI = sprintf('%s_SI_%d_%.2f_00000%d.bmp', SeqName, GOPSize, MKP, fn - 1);
        elseif (fn - 1 >= 10) && (fn - 1 <= 99)
            imgfile = sprintf('%s_0000%d.bmp', SeqName, fn - 1);
            Recimgfile = sprintf('Rec%s0%d.bmp', SeqName, fn - 1); 
            SI = sprintf('%s_SI_%d_%.2f_0000%d.bmp', SeqName, GOPSize, MKP, fn - 1);
        else
            imgfile = sprintf('%s_000%d.bmp', SeqName, fn - 1);
            Recimgfile = sprintf('Rec%s%d.bmp', SeqName, fn - 1);
            SI = sprintf('%s_SI_%d_%.2f_000%d.bmp', SeqName, GOPSize, MKP, fn - 1);
        end                
        
        [Rec_CS, CS_PSNR, num_nz_theta, Num_Iteration] = My_fast_cs2d(imgfile, MCS, 'BWHT', 0, 32, SI, Iteration);
        Ave_num_nz_theta = Ave_num_nz_theta + num_nz_theta;
        Ave_PSNR = Ave_PSNR + CS_PSNR;
        Ave_CS_PSNR = Ave_CS_PSNR + CS_PSNR;
        Ave_num_iteration = Ave_num_iteration + Num_Iteration;
        imwrite(Rec_CS, Recimgfile, 'bmp');
        CSFN = CSFN + 1;
    end
end

Total_time = toc; % total time
fprintf(1, 'Number of reconstructed KEY frames is %d\n', KFN);
fprintf(1, 'Number of reconstructed CS frames is %d\n', CSFN);
Ave_Obr_Rate = (KFN * MK + CSFN * MCS) / (NF * FrameWidth * FrameHeight);
fprintf(1, 'Average measurement rate per frame is %.2f\n', Ave_Obr_Rate);
Ave_K_PSNR = Ave_K_PSNR / KFN;
fprintf(1, 'Average PSNR of the reconstructed KEY frames is %.2fdB\n', Ave_K_PSNR);
Ave_CS_PSNR = Ave_CS_PSNR / CSFN;
fprintf(1, 'Average PSNR of the reconstructed CS frames is %.2fdB\n', Ave_CS_PSNR);
Ave_PSNR = Ave_PSNR / NF;
fprintf(1, 'Average PSNR of the reconstructed video sequence is %.2fdB\n', Ave_PSNR);
Ave_num_nz_theta = Ave_num_nz_theta / NF;
fprintf(1, 'Average number of non-zero components of the reconstructed video sequence is %.2f\n', Ave_num_nz_theta);
fprintf(1, 'Average rate of non-zero components of the reconstructed video sequence is %.4f\n', Ave_num_nz_theta / (FrameWidth * FrameHeight));
Ave_num_iteration = Ave_num_iteration / NF;
fprintf(1, 'Average number of iteration per frame is %.2f iterations\n', Ave_num_iteration);
fprintf(1, 'Average reconstruction time per frame is %.2f seconds\n\n', (Total_time + 10) / NF);
close all;