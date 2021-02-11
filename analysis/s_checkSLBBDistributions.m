bins = 20;

ntest = @(x) kstest((x-mean(x, 'omitnan'))/std(x,'omitnan'));

sl1 = ampl{s}.sl.full(:,14);
bb1 = sqrt(ampl{s}.bb.full(:,14));
% bb1_blank = sqrt(ampl{s}.bb.blank(:,2));

[h,p]=ntest(sl1)
[h,p]=ntest(log(sl1))
[h,p]=ntest(sl1.^2)
[h,p]=ntest(log(sl1.^2))

[h,p]=ntest(bb1)
[h,p]=ntest(log(bb1))
[h,p]=ntest(bb1.^2)
[h,p]=ntest(log(bb1.^2))

figure;
subplot(2,2,1);
histogram(sl1, bins)
title('amplitudes')

subplot(2,2,2);
histogram(log(sl1), bins)
title('log amplitudes')

subplot(2,2,3);
histogram(sl1.^2, bins)
title('power ')

subplot(2,2,4);
histogram(log(sl1.^2), bins)
title('log power')


figure;
subplot(2,2,1);
histogram(bb1, bins)
title('bb amplitudes')

subplot(2,2,2);
histogram(log(bb1), bins)
title('bb log amplitudes')

subplot(2,2,3);
histogram(bb1.^2, bins)
title('bb power ')

subplot(2,2,4);
histogram(log(bb1.^2), bins)
title('bb log power')

