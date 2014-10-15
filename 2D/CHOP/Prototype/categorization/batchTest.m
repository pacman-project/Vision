accuracyArr = [6,1];
for itr = 3:6
    TrainCategoryModel('CIFAR10', 1, itr);
    accuracy = TestCategoryModel('CIFAR10', 1, itr);
    accuracyArr(itr) = accuracy;
end