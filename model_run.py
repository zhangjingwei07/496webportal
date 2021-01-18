class model_run:
    ###------------------------------Model---------------------------------------###
    TOP_K = 5
    # 计算Binary和得到labels
    def compute_result(dataloader, net):
        binariy_codes, labels = [], []
        net.eval()
        for img, label in dataloader:
            labels.append(label)
            binariy_codes.append((net(img.cuda())).data)
        return torch.cat(binariy_codes).sign(), torch.cat(labels)


    class TMDataModule(pl.LightningDataModule):
      def __init__(self, data_dir: str = './data', batch_size=64, num_workers=2, import_size=0.4):
        super().__init__()
        self.data_dir = data_dir
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.data_name = "TM"
        self.import_size = import_size # Percentage of original dataset, Total dataset size = 54865. Went to 0.4 without crashing
        self.label_mapping = None
        self.N_FEATURES = 19791
        self.N_CLASS = 55

      def prepare_data(self):
        # download
        if not os.path.exists(f"{self.data_dir}/{self.data_name}"):
          print("Start downloading TM!")
          url = 'https://github.com/Aprilhuu/Deep-Learning-in-Single-Cell-Analysis/raw/main/TM.zip'
          download_and_extract_archive(url, f"{self.data_dir}", f"{self.data_dir}", filename=self.data_name + '.zip')

      def setup(self, stage):
        # Assign train/val datasets for use in dataloaders

        DataPath = self.data_dir + "/" + self.data_name + "/Filtered_TM_data.csv"
        LabelsPath = self.data_dir + "/" + self.data_name + "/Labels.csv"
        
        # Step #1: Read in all labels
        labels = pd.read_csv(LabelsPath, header=0, index_col=None, sep=',')
        full_labels = np.asarray(labels)
        labels = None

        # Step #2: Preprocess data and only keep cells with population larger than 10
        remaining_labels, discarded_indices = preprocess_data(full_labels, self.import_size)

        # Step #3: Turning string class labels into int labels and store the mapping
        self.label_mapping = preprocessing.LabelEncoder()
        self.label_mapping.fit(remaining_labels)
        int_labels = self.label_mapping.transform(remaining_labels)

        full_labels = np.asarray(int_labels)
        remaining_labels = None
        int_labels = None

        # Step #4: Read in data based on selected label indices
        discarded_indices = [x + 1 for x in discarded_indices]
        data = pd.read_csv(DataPath, index_col=0, sep=',', skiprows=discarded_indices)
        discarded_indices = None
        full_data = np.asarray(data, dtype=np.float32)
        full_dataset = CustomDataset(full_data, full_labels)

        # Step #5: Split indices in stratified way into train, validation, test and database sets
        #          and prepare all datasets based on splited list indices
        self.TM_database, self.TM_train, self.TM_val, self.TM_test = split_test_train_val_database_sets(full_dataset,
                                                                                                        train_percentage=0.3,
                                                                                                        val_percentage=0.1,
                                                                                                        test_percentage=0.1)
        
        self.database_dataloader = DataLoader(self.TM_database, batch_size=self.batch_size,
                                              num_workers=self.num_workers)
         
        print("database size =", len(self.TM_database))
        print("train size =", len(self.TM_train))
        print("val size =", len(self.TM_val))
        print("test size =", len(self.TM_test))
        
        # Calculate sample count in each class for training dataset
        samples_in_each_class_dict = Counter([data[1] for data in self.TM_train])
        print("training samples in each class =", sorted(samples_in_each_class_dict.items()))
        samples_in_each_class_dict_val = Counter([data[1] for data in self.TM_val])
        print("val samples in each class =", sorted(samples_in_each_class_dict_val.items()))
        samples_in_each_class_dict_test = Counter([data[1] for data in self.TM_test])
        print("test samples in each class =", sorted(samples_in_each_class_dict_test.items()))
        self.N_CLASS = len(samples_in_each_class_dict)
        print("Changing N_CLASS =", self.N_CLASS)
        self.samples_in_each_class = torch.zeros(self.N_CLASS)
        for index, count in samples_in_each_class_dict.items():
           self.samples_in_each_class[index] = count


    def train_dataloader(self):
        return DataLoader(self.TM_train, batch_size=self.batch_size,
                          shuffle=True, num_workers=self.num_workers)

    def val_dataloader(self):
        return DataLoader(self.TM_val, batch_size=self.batch_size,
                          num_workers=self.num_workers)

    def test_dataloader(self):
        return DataLoader(self.TM_test, batch_size=self.batch_size,
                          num_workers=self.num_workers)
    class CSQLightening(pl.LightningModule):
      def __init__(self,n_class,n_features,batch_size=64,l_r=1e-5,lamb_da=0.0001,beta=0.9999,bit=64):
        super(CSQLightening, self).__init__()
        print("hparam: l_r = {}, lambda = {}, beta = {}", l_r, lamb_da, beta)
        self.batch_size = batch_size
        self.l_r = l_r
        self.bit = bit
        self.n_class = n_class
        self.lamb_da = lamb_da
        self.beta = beta
        self.samples_in_each_class = None # Later initialized in training step
        self.hash_centers = get_hash_centers(self.n_class, self.bit)
        ##### model structure ####
        # input size = batch size * 19791
        self.hash_layer = nn.Sequential(
            nn.Linear(n_features, 6300),
            nn.ReLU(inplace=True),
            nn.Dropout(0.2),
            nn.Linear(6300, 2100),
            nn.ReLU(inplace=True),
            nn.Dropout(0.2),
            nn.Linear(2100, 710),
            nn.ReLU(inplace=True),
            nn.Dropout(0.2),
            nn.Linear(710, 200),
            nn.ReLU(inplace=True),
            nn.Linear(200, self.bit),
        )
        

      def forward(self, x):
        # forward pass returns prediction
          x = self.hash_layer(x)
          return x

      def CSQ_loss_function(self, hash_codes, labels):
        hash_codes = hash_codes.tanh()
        hash_centers = self.hash_centers[labels]
        hash_centers = hash_centers.type_as(hash_codes)

        # Class-Balanced Loss on Effective Number of Samples
        # Reference Paper https://arxiv.org/abs/1901.05555
        if self.samples_in_each_class == None:
          self.samples_in_each_class = self.trainer.datamodule.samples_in_each_class
          self.n_class = self.trainer.datamodule.N_CLASS
        class_sample_count = self.samples_in_each_class[labels]
        weight = (1 - self.beta)/(1 - torch.pow(self.beta, class_sample_count))
        weight = weight / weight.sum() * self.n_class
        weight = weight.type_as(hash_codes)

        # Center Similarity Loss
        BCELoss = nn.BCELoss(weight=weight.unsqueeze(1).repeat(1,self.bit))
        C_loss = BCELoss(0.5 * (hash_codes + 1),
                            0.5 * (hash_centers + 1))
        # Quantization Loss
        Q_loss = (hash_codes.abs() - 1).pow(2).mean()

        loss = C_loss + self.lamb_da * Q_loss
        return loss


    ## Initialize existing model
    def load_model(PATH):
        model = CSQLightening.load_from_checkpoint(checkpoint_path=PATH)
        
        return model

    # Predict label using Closest Hash Center strategy (b)
    def get_labels_pred_closest_hash_center(query_binaries, query_labels, hash_centers):
        num_query = query_labels.shape[0]
        labels_pred = []
        for binary_query, label_query in zip(query_binaries, query_labels):
              dists = CalcHammingDist(binary_query, hash_centers)
              closest_class = np.argmin(dists)
              labels_pred.append(closest_class)
        return labels_pred

    def run(query_dataloader):
        model = load_model('./checkpoint.ckpt')
        datamodule = TMDataModule(import_size=0.1, num_workers=4)
        
        binaries_database, labels_database = compute_result(model.trainer.datamodule.database_dataloader, net)
        binaries_query, labels_query = compute_result(query_dataloader, model)

        # 转换成one-hot encoding，方便后续高效计算
        labels_database_one_hot = categorical_to_onehot(labels_database, class_num)
        labels_query_one_hot = categorical_to_onehot(labels_query, class_num)
        labels_pred_CHC = get_labels_pred_closest_hash_center(binaries_query.cpu().numpy(), labels_query.numpy(), net.hash_centers.numpy())
        
        return labels_pred_CHC
