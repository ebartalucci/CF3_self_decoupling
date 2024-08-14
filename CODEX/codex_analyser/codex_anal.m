FOLDER = rbnmr('D:\PhD\Data\NMR\ETBA\Bolm\CF3_coupling\130824_3p2mm_19F_TLA_Static_and_CODEX');
ExpNoStart = 200;
ExpNoEnd = 231;
NoPeaks = 1;
plotbnmr(FOLDER(13:44))








% Function section
function RMSD = computeRMSD(I_sim, I_exp, N)
    RMSD = sqrt(sum((I_sim - I_exp)/N));
end

function k_ii = computeRateConstant(r_ij, theta_ij)

end

