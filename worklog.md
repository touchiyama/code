#
# RaptGen Worklog
#

## * ToDo
* For embed size, change 2(default) to 10, which leads to get the center of the GMM populations via gmm.py.
* Using GPU on Google colaboratory, conduct decoder.py.


## * SetUp
```bash
TomoyaUchiyama@MACOSMBP1902 RaptGen % brew install cairo && brew install pango
```
```bash
TomoyaUchiyama@MACOSMBP1902 RaptGen % pip install biopython 
TomoyaUchiyama@MACOSMBP1902 RaptGen % pip install CairoSVG svgwrite
TomoyaUchiyama@MACOSMBP1902 RaptGen % pip install GPyOpt==1.2.6
```

```bash
TomoyaUchiyama@MACOSMBP1902 RaptGen % conda create -n RaptGen python=3.8
TomoyaUchiyama@MACOSMBP1902 RaptGen % pyenv local anaconda3-2023.03/envs/RaptGen
TomoyaUchiyama@MACOSMBP1902 RaptGen % pyenv versions                            
  system
  anaconda3-2022.05
  anaconda3-2023.03
* anaconda3-2023.03/envs/RaptGen (set by /Users/TomoyaUchiyama/RaptGen/.python-version)
```
```bash
TomoyaUchiyama@MACOSMBP1902 RaptGen % python3 setup.py install
TomoyaUchiyama@MACOSMBP1902 RaptGen % pip install click==7.1.2
TomoyaUchiyama@MACOSMBP1902 RaptGen % pip install numpy
TomoyaUchiyama@MACOSMBP1902 RaptGen % pip install torch==1.5.0
TomoyaUchiyama@MACOSMBP1902 RaptGen % pip install tqdm==4.41.1
TomoyaUchiyama@MACOSMBP1902 RaptGen % pip install matplotlib==3.1.2
TomoyaUchiyama@MACOSMBP1902 RaptGen % pip install pandas
TomoyaUchiyama@MACOSMBP1902 RaptGen % pip install Pillow
TomoyaUchiyama@MACOSMBP1902 RaptGen % pip install scikit-learn
TomoyaUchiyama@MACOSMBP1902 RaptGen % pip install Pathlib
```
```bash
# check script ---
TomoyaUchiyama@MACOSMBP1902 RaptGen % more setup.py
from setuptools import find_packages, setup

setup(
    name='raptgen',
    packages=find_packages(),
    version='0.1.0',
    description='A short description of the project.',
    author='Natsuki Iwano',
    license='MIT',
)
```
## * Make sure to run rightly
### (0) compile scripts modified for aa analysis.
```bash
TomoyaUchiyama@MACOSMBP1902 RaptGen % python3 setup.py install

running install
/Users/TomoyaUchiyama/.pyenv/versions/anaconda3-2023.03/envs/RaptGen/lib/python3.8/site-packages/setuptools/_distutils/cmd.py:66: SetuptoolsDeprecationWarning: setup.py install is deprecated.
!!

        ********************************************************************************
        Please avoid running ``setup.py`` directly.
        Instead, use pypa/build, pypa/installer, pypa/build or
        other standards-based tools.

        See https://blog.ganssle.io/articles/2021/10/setup-py-deprecated.html for details.
        ********************************************************************************

!!
  self.initialize_options()
/Users/TomoyaUchiyama/.pyenv/versions/anaconda3-2023.03/envs/RaptGen/lib/python3.8/site-packages/setuptools/_distutils/cmd.py:66: EasyInstallDeprecationWarning: easy_install command is deprecated.
!!

        ********************************************************************************
        Please avoid running ``setup.py`` and ``easy_install``.
        Instead, use pypa/build, pypa/installer, pypa/build or
        other standards-based tools.

        See https://github.com/pypa/setuptools/issues/917 for details.
        ********************************************************************************

!!
  self.initialize_options()
running bdist_egg
running egg_info
writing raptgen.egg-info/PKG-INFO
writing dependency_links to raptgen.egg-info/dependency_links.txt
writing top-level names to raptgen.egg-info/top_level.txt
reading manifest file 'raptgen.egg-info/SOURCES.txt'
adding license file 'LICENSE'
writing manifest file 'raptgen.egg-info/SOURCES.txt'
installing library code to build/bdist.macosx-11.1-arm64/egg
running install_lib
running build_py
copying raptgen/visualization_aa.py -> build/lib/raptgen
copying raptgen/models_aa.py -> build/lib/raptgen
copying raptgen/data.py -> build/lib/raptgen
copying raptgen/data_aa.py -> build/lib/raptgen
creating build/bdist.macosx-11.1-arm64/egg
creating build/bdist.macosx-11.1-arm64/egg/raptgen
copying build/lib/raptgen/models.py -> build/bdist.macosx-11.1-arm64/egg/raptgen
copying build/lib/raptgen/visualization_aa.py -> build/bdist.macosx-11.1-arm64/egg/raptgen
copying build/lib/raptgen/models_aa.py -> build/bdist.macosx-11.1-arm64/egg/raptgen
copying build/lib/raptgen/__init__.py -> build/bdist.macosx-11.1-arm64/egg/raptgen
copying build/lib/raptgen/visualization.py -> build/bdist.macosx-11.1-arm64/egg/raptgen
copying build/lib/raptgen/data.py -> build/bdist.macosx-11.1-arm64/egg/raptgen
copying build/lib/raptgen/data_aa.py -> build/bdist.macosx-11.1-arm64/egg/raptgen
byte-compiling build/bdist.macosx-11.1-arm64/egg/raptgen/models.py to models.cpython-38.pyc
byte-compiling build/bdist.macosx-11.1-arm64/egg/raptgen/visualization_aa.py to visualization_aa.cpython-38.pyc
byte-compiling build/bdist.macosx-11.1-arm64/egg/raptgen/models_aa.py to models_aa.cpython-38.pyc
byte-compiling build/bdist.macosx-11.1-arm64/egg/raptgen/__init__.py to __init__.cpython-38.pyc
byte-compiling build/bdist.macosx-11.1-arm64/egg/raptgen/visualization.py to visualization.cpython-38.pyc
byte-compiling build/bdist.macosx-11.1-arm64/egg/raptgen/data.py to data.cpython-38.pyc
byte-compiling build/bdist.macosx-11.1-arm64/egg/raptgen/data_aa.py to data_aa.cpython-38.pyc
creating build/bdist.macosx-11.1-arm64/egg/EGG-INFO
copying raptgen.egg-info/PKG-INFO -> build/bdist.macosx-11.1-arm64/egg/EGG-INFO
copying raptgen.egg-info/SOURCES.txt -> build/bdist.macosx-11.1-arm64/egg/EGG-INFO
copying raptgen.egg-info/dependency_links.txt -> build/bdist.macosx-11.1-arm64/egg/EGG-INFO
copying raptgen.egg-info/top_level.txt -> build/bdist.macosx-11.1-arm64/egg/EGG-INFO
zip_safe flag not set; analyzing archive contents...
creating 'dist/raptgen-0.1.0-py3.8.egg' and adding 'build/bdist.macosx-11.1-arm64/egg' to it
removing 'build/bdist.macosx-11.1-arm64/egg' (and everything under it)
Processing raptgen-0.1.0-py3.8.egg
Removing /Users/TomoyaUchiyama/.pyenv/versions/anaconda3-2023.03/envs/RaptGen/lib/python3.8/site-packages/raptgen-0.1.0-py3.8.egg
Copying raptgen-0.1.0-py3.8.egg to /Users/TomoyaUchiyama/.pyenv/versions/anaconda3-2023.03/envs/RaptGen/lib/python3.8/site-packages
Adding raptgen 0.1.0 to easy-install.pth file

Installed /Users/TomoyaUchiyama/.pyenv/versions/anaconda3-2023.03/envs/RaptGen/lib/python3.8/site-packages/raptgen-0.1.0-py3.8.egg
Processing dependencies for raptgen==0.1.0
Finished processing dependencies for raptgen==0.1.0
```
### (1) Build up model (training a.a seq generation model) 
```bash
TomoyaUchiyama@MACOSMBP1902 RaptGen % python3 /Users/TomoyaUchiyama/RaptGen/scripts/real_aa.py --help 
Usage: real_aa.py [OPTIONS] SEQPATH

  run experiment with real data

Options:
  --epochs INTEGER        the number of training epochs  [default: 1000]
  --threshold INTEGER     the number of epochs with no loss update to stop
                          training  [default: 50]

  --use-cuda / --no-cuda  use cuda if available  [default: True]
  --cuda-id INTEGER       the device id of cuda to run  [default: 0]
  --save-dir PATH         path to save results  [default:
                          /Users/TomoyaUchiyama/RaptGen/out/real]

  --fwd TEXT              forward adapter
  --rev TEXT              reverse adapter
  --min-count INTEGER     minimum duplication count to pass sequence for
                          training  [default: 1]

  --multi INTEGER         the number of training for multiple times  [default:
                          1]

  --reg-epochs INTEGER    the number of epochs to conduct state transition
                          regularization  [default: 50]

  --embed-size INTEGER    the number of embedding dimension of raptgen model
                          [default: 2]

  --fast / --normal       [experimental] use fast calculation of probability
                          estimation. Output of the decoder shape is different
                          and the visualizers are not implemented.  [default:
                          False]

  --help                  Show this message and exit.  [default: False]
```
```bash
TomoyaUchiyama@MACOSMBP1902 RaptGen % python3 /Users/TomoyaUchiyama/RaptGen/scripts/real_aa.py \
/Users/TomoyaUchiyama/RaptGen/data/simulation/paired/sequences_aa.fasta \
--save-dir /Users/TomoyaUchiyama/RaptGen/out_aa

saving to /Users/TomoyaUchiyama/RaptGen/out_aa
reading fasta format sequence
adapter info not provided. estimating value
no match found.
no match found.
no match found.
no match found.
filtering with : (0N)-6N-(0N)
filtering with : (0N)-6N-(0N)
experiment name : 20230702_201251240308
# of sequences -> 4187
# of sequences -> 4187
training cnn_phmm_vae.mdl
[50] 234 itr  11.31 <->  11.34 ( 10.55+  0.79) of .cnn_phmm_:  23%|█████████▎                              | 234/1000 [41:16<2:15:07, 10.58s/it]
```
* embed size is 10 
```bash
(base) TomoyaUchiyama@MACOSMBP1902 RaptGen % python3 /Users/TomoyaUchiyama/RaptGen/scripts/real_aa.py \
/Users/TomoyaUchiyama/RaptGen/data/simulation/paired/sequences_aa.fasta \
--save-dir /Users/TomoyaUchiyama/RaptGen/out_aa \
--embed-size 10

(base) TomoyaUchiyama@MACOSMBP1902 RaptGen % cd out_aa
(base) TomoyaUchiyama@MACOSMBP1902 out_aa % mv cnn_phmm_vae.csv cnn_phmm_vae_10.csv
(base) TomoyaUchiyama@MACOSMBP1902 out_aa % mv cnn_phmm_vae.mdl cnn_phmm_vae_10.mdl
```

### (2) Encoder (sequence vector in embedded space)
```bash
TomoyaUchiyama@MACOSMBP1902 RaptGen % python3 /Users/TomoyaUchiyama/RaptGen/scripts/encode_aa.py --help
Usage: encode_aa.py [OPTIONS] SEQPATH MODELPATH

  achieve sequence vector in embedded space.

Options:
  --use-cuda / --no-cuda  use cuda if available  [default: True]
  --cuda-id INTEGER       the device id of cuda to run  [default: 0]
  --fwd TEXT              forward adapter
  --rev TEXT              reverse adapter
  --save-dir PATH         path to save results  [default:
                          /Users/TomoyaUchiyama/RaptGen/out/encode]

  --fast / --normal       [experimental] use fast calculation of probability
                          estimation.         Output of the decoder shape is
                          different and the visualizers are not implemented.
                          [default: False]

  --help                  Show this message and exit.  [default: False]
```
```bash
TomoyaUchiyama@MACOSMBP1902 RaptGen % python3 /Users/TomoyaUchiyama/RaptGen/scripts/encode_aa.py \
/Users/TomoyaUchiyama/RaptGen/data/simulation/paired/sequences_aa.fasta \
/Users/TomoyaUchiyama/RaptGen/out_aa/cnn_phmm_vae.mdl \
--save-dir /Users/TomoyaUchiyama/RaptGen/out_aa
saving to /Users/TomoyaUchiyama/RaptGen/out_aa
reading fasta format sequence
adapter info not provided. estimating value
no match found.
no match found.
no match found.
no match found.
filtering with : (0N)-6N-(0N)
filtering with : (0N)-6N-(0N)
experiment name : 20230702_210403048476
skip loading training result
skip loading model parameters
evaluating mu
# of sequences -> 4187
# of sequences -> 4187
saving to /Users/TomoyaUchiyama/RaptGen/out_aa/embed_seq.csv
... done.
```bash


```


```
### (3) Decoder (latent point to most_probable sequence)
```bash
TomoyaUchiyama@MACOSMBP1902 RaptGen % python3 /Users/TomoyaUchiyama/RaptGen/scripts/decode_aa.py --help
Usage: decode_aa.py [OPTIONS] POS_PATH MODEL_PATH TARGET_LEN

  achieve sequence vector in embedded space.

Options:
  --use-cuda / --no-cuda  use cuda if available  [default: True]
  --cuda-id INTEGER       the device id of cuda to run  [default: 0]
  --save-dir PATH         path to save results  [default:
                          /Users/TomoyaUchiyama/RaptGen/out/decode]

  --embed-dim INTEGER     the embedding dimension of raptgen model  [default:
                          2]

  --eval-max INTEGER      the maximum number of sequence to evaluate most
                          probable sequence  [default: 256]

  --help                  Show this message and exit.  [default: False]
```
```bash
TomoyaUchiyama@MACOSMBP1902 RaptGen % python3 /Users/TomoyaUchiyama/RaptGen/scripts/decode_aa.py \
/Users/TomoyaUchiyama/RaptGen/out_aa/embed_seq.csv \
/Users/TomoyaUchiyama/RaptGen/out_aa/cnn_phmm_vae.mdl \
6
--save-dir /Users/TomoyaUchiyama/RaptGen/out_aa
saving to /Users/TomoyaUchiyama/RaptGen/out/decode
loading model parameters
calculating phmm parameter
generating sequences
0: ('*TT*', 'ATTH', -4.485268592834473)
1: ('*LL*', 'FLLM', -6.542054176330566)
2: ('*PS*', 'YPSM', -6.533905506134033)
3: ('*TL*', 'HTLC', -5.667392253875732)
4: ('*PS*', 'KPSI', -5.679365634918213)
```

### (4) gmm
* select the center of the GMM populations
* output the top 10 sequences to a specified directory
```bash
# gmm.py ---
(base) TomoyaUchiyama@MACOSMBP1902 RaptGen % python3 /Users/TomoyaUchiyama/RaptGen/scripts/gmm_aa.py --help
Usage: gmm_aa.py [OPTIONS] SEQPATH MODELPATH

  select gmm center with trained model

Options:
  --use-cuda / --no-cuda  use cuda if available  [default: use-cuda]
  --cuda-id INTEGER       the device id of cuda to run  [default: 0]
  --save-dir PATH         path to save results  [default:
                          /Users/TomoyaUchiyama/RaptGen/out/gmm]
  --fwd TEXT              forward adapter
  --rev TEXT              reverse adapter
  --help                  Show this message and exit.  [default: False]
```
```bash
(base) TomoyaUchiyama@MACOSMBP1902 RaptGen % python3 /Users/TomoyaUchiyama/RaptGen/scripts/gmm_aa.py \
/Users/TomoyaUchiyama/RaptGen/data/simulation/paired/sequences_aa.fasta \
/Users/TomoyaUchiyama/RaptGen/out_aa/cnn_phmm_vae_10.mdl \
--save-dir /Users/TomoyaUchiyama/RaptGen/out_aa
```

## * Note
* fonts-noto-cjk
```bash
Noto CJK fonts are rebranded versions of Adobe Source Han fonts, developed by Adobe and Google which contains Chinese characters, Hangul and Kana; Latin-script letters and numerals are taken from the Source Pro fonts.
```

* Train model 
```bash
from raptgen.models import CNN_PHMM_VAE, CNN_PHMM_VAE_FAST
from raptgen.data import SequenceGenerator, SingleRound
# -> remove SequenceGenerator
```
* CNN_PHMM_VAE <br>
-> VAE is an inherited class.
```bash
class CNN_PHMM_VAE(VAE):
    def __init__(self, motif_len=12, embed_size=10, hidden_size=32, kernel_size=7):
        encoder = EncoderCNN(hidden_size, kernel_size)
        decoder = DecoderPHMM(motif_len, embed_size)

        super(CNN_PHMM_VAE, self).__init__(
            encoder, decoder, embed_size, hidden_size)
        self.loss_fn = profile_hmm_loss_fn
```
```bash
class VAE(nn.Module):
    def __init__(self, encoder, decoder, embed_size=10, hidden_size=32):
        super(VAE, self).__init__()

        self.encoder = encoder
        self.decoder = decoder

        self.h2mu = nn.Linear(hidden_size, embed_size)
        self.h2logvar = nn.Linear(hidden_size, embed_size)

    def reparameterize(self, mu, logvar, deterministic=False):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        z = mu + (std * eps if not deterministic else 0)
        return z

    def forward(self, input, deterministic=False):
        h = self.encoder(input)
        mu = self.h2mu(h)
        logvar = self.h2logvar(h)

        z = self.reparameterize(mu, logvar, deterministic)
        recon_param = self.decoder(z)
        return recon_param, mu, logvar
```
```bash
class Bottleneck(nn.Module):
    def __init__(self, init_dim=32, window_size=7):
        super(Bottleneck, self).__init__()
        assert window_size % 2 == 1, f"window size should be odd, given {window_size}"

        self.conv1 = nn.Conv1d(
            in_channels=init_dim,
            out_channels=init_dim*2,
            kernel_size=1)

        self.conv2 = nn.Conv1d(
            in_channels=init_dim*2,
            out_channels=init_dim*2,
            kernel_size=window_size,
            padding=window_size//2
        )

        self.conv3 = nn.Conv1d(
            in_channels=init_dim*2,
            out_channels=init_dim,
            kernel_size=1)

        self.bn1 = nn.BatchNorm1d(init_dim)
        self.bn2 = nn.BatchNorm1d(init_dim*2)
        self.bn3 = nn.BatchNorm1d(init_dim*2)

    def forward(self, input):
        x = self.conv1(F.leaky_relu(self.bn1(input)))
        x = self.conv2(F.leaky_relu(self.bn2(x)))
        x = self.conv3(F.leaky_relu(self.bn3(x)))
        return F.leaky_relu(x+input)
```
```bash
class EncoderCNN (nn.Module):
    # 0~3 is already used by embedding ATGC
    def __init__(self, embedding_dim=32, window_size=7, num_layers=6):
        super(EncoderCNN, self).__init__()
        self.embedding_dim = embedding_dim
        self.window_size = window_size

        self.embed = nn.Embedding(
            num_embeddings=4,  # [A,T,G,C,PAD,SOS,EOS]
            embedding_dim=embedding_dim)

        modules = [Bottleneck(embedding_dim, window_size)
                   for _ in range(num_layers)]
        self.resnet = nn.Sequential(*modules)

    def forward(self, seqences):
        # change X from (N, L) to (N, L, C)
        x = F.leaky_relu(self.embed(seqences))

        # change X to (N, C, L)
        x = x.transpose(1, 2)
        value, indices = self.resnet(x).max(dim=2)
        return value
```
```bash
class DecoderPHMM(nn.Module):
    # tile hidden and input to make x
    def __init__(self,  motif_len, embed_size,  hidden_size=32):
        super(DecoderPHMM, self).__init__()

        class View(nn.Module):
            def __init__(self, shape):
                super(View, self).__init__()
                self.shape = shape

            def forward(self, x):
                return x.view(*self.shape)

        self.fc1 = nn.Sequential(
            nn.Linear(embed_size, hidden_size),
            nn.BatchNorm1d(hidden_size),
            nn.LeakyReLU(negative_slope=0.01, inplace=True))

        self.fc2 = nn.Sequential(
            nn.Linear(hidden_size, hidden_size*2),
            nn.BatchNorm1d(hidden_size*2),
            nn.LeakyReLU(negative_slope=0.01),
            nn.Linear(hidden_size*2, hidden_size),
            nn.BatchNorm1d(hidden_size),
            nn.LeakyReLU(negative_slope=0.01)
        )

        self.tr_from_M = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.LeakyReLU(negative_slope=0.01, inplace=True),
            nn.Linear(hidden_size, (motif_len+1)*3),
            nn.LeakyReLU(negative_slope=0.01, inplace=True),
            View((-1, motif_len+1, 3)),
            nn.LogSoftmax(dim=2)
        )
        self.tr_from_I = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.LeakyReLU(negative_slope=0.01, inplace=True),
            nn.Linear(hidden_size, (motif_len+1)*2),
            nn.LeakyReLU(negative_slope=0.01, inplace=True),
            View((-1, motif_len+1, 2)),
            nn.LogSoftmax(dim=2)
        )
        self.tr_from_D = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.LeakyReLU(negative_slope=0.01, inplace=True),
            nn.Linear(hidden_size, (motif_len+1)*2),
            nn.LeakyReLU(negative_slope=0.01, inplace=True),
            View((-1, motif_len+1, 2)),
            nn.LogSoftmax(dim=2)
        )

        self.emission = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.LeakyReLU(negative_slope=0.01, inplace=True),
            nn.Linear(hidden_size, motif_len*4),
            nn.LeakyReLU(negative_slope=0.01, inplace=True),
            View((-1, motif_len, 4)),
            nn.LogSoftmax(dim=2)
        )

    def forward(self, input):
        x = self.fc1(input)

        transition_from_match = self.tr_from_M(x)
        transition_from_insertion = self.tr_from_I(x)
        transition_from_deletion = self.tr_from_D(x)

        emission_proba = self.emission(x)
        return (torch.cat((
            transition_from_match,
            transition_from_insertion,
            transition_from_deletion), dim=2), emission_proba)
```
```bash
class SingleRound:
    """pass path or raw_reads to make class of selex experiment per round.
    """

    def __init__(self, raw_reads: list = None, forward_adapter=None, reverse_adapter=None, name=None, tolerance=0, path: str = None, max_len=None):
        
        assert path is not None or raw_reads is not None, \
        "either path or raw_reads has to be specified"

        if path:
            path = Path(path)
            if path .suffix == ".fastq":
                logger.info("reading fastq format sequence")
                raw_reads = read_fastq(path)
            elif path.suffix in {".fasta", ".fa"}:
                logger.info("reading fasta format sequence")
                raw_reads = read_fasta(path)
            else:
                logger.critical(
                    "please specify a file with fasta or fastq format")
                quit()

        self.raw_reads = raw_reads
        self.calc_target_length()
        self.max_len = max_len

        if forward_adapter is None or reverse_adapter is None:
            logger.info("adapter info not provided. estimating value")
            self.calc_experimental_settings()
        
        else:
            logger.info(
                f"sequence design : \
                {forward_adapter}-[random]-{reverse_adapter}"
            )
            self.set_adapters(forward_adapter, reverse_adapter,
                              self.max_len is not None)

        if name:
            self.name = name
        else:
            self.name = re.sub(
                r'[-\.\:]', 
                "",
                str(datetime.datetime.now())
            ).replace(" ", "_")
        
        logger.info(f"experiment name : {self.name}")
        self.tolerance = tolerance

    def get_adapters(self):
        return self.forward_adapter, self.reverse_adapter

    def set_adapters(self, forward_adapter: str, reverse_adapter: str, set_max_len=False):
        self.forward_adapter = forward_adapter
        self.forward_adapter_length = len(forward_adapter)

        self.reverse_adapter = reverse_adapter
        self.reverse_adapter_length = len(reverse_adapter)

        self.random_region_length = self.target_length - \
            self.reverse_adapter_length - self.forward_adapter_length
        if set_max_len:
            self.random_region_length = self.max_len

    def calc_target_length(self):
        from collections import Counter, defaultdict
        self.read_counter = Counter(self.raw_reads)

        # calc most common length
        d = defaultdict(int)
        for key, value in self.read_counter.items():
            d[len(key)] += value
        self.target_length = sorted(d.items(), key=lambda x: -x[1])[0][0]

    def calc_experimental_settings(self):
        """calculate sequence adapters in a heuristic way
        """

        # fwd
        max_count = None
        est_adapter = ""
        for i in range(1, self.target_length):
            d = defaultdict(int)
            for seq, count in self.read_counter.most_common():
                if len(seq) < i or len(d) > 100 and seq[:i] not in d.keys():
                    continue
                d[seq[:i]] += count
            top_seq, top_count = sorted(d.items(), key=lambda x: -x[1])[0]
            if max_count is not None and top_count < max_count * 0.5:  # heuristics
                logger.info(
                    f"estimated forward adapter len is {i-1} : {est_adapter}")
                break
            max_count = sorted(d.items(), key=lambda x: -x[1])[0][1]
            if max_count < sum(self.read_counter.values()) * 0.5:
                logger.info(
                    f"no match found.")
                break
            est_adapter = top_seq
        fwd_len = i - 1
        fwd_adapter = est_adapter

        # rev
        max_count = None
        est_adapter = ""
        for i in range(1, self.target_length):
            d = defaultdict(int)
            for seq, count in self.read_counter.most_common():
                if len(seq) < i or len(d) > 100 and seq[-i:] not in d.keys():
                    continue
                d[seq[-i:]] += count
            top_seq, top_count = sorted(d.items(), key=lambda x: -x[1])[0]
            if max_count is not None and top_count < max_count * 0.5:  # heuristics
                logger.info(
                    f"estimated reverse adapter len is {i-1} : {est_adapter}")
                break
            max_count = sorted(d.items(), key=lambda x: -x[1])[0][1]
            if max_count < sum(self.read_counter.values()) * 0.5:
                logger.info(
                    f"no match found.")
                break
            est_adapter = top_seq
        rev_len = i - 1
        rev_adapter = est_adapter

        rand_len = self.target_length - rev_len - fwd_len

        logger.info(
            f"filtering with : {fwd_adapter}({fwd_len}N)-{rand_len}N-{rev_adapter}({rev_len}N)")

        # write estimated experimental settings
        self.set_adapters(fwd_adapter, rev_adapter, self.max_len is not None)

    def get_sequences_and_count(self):
        c = Counter(self.raw_reads)
        return c.most_common()

    def get_filter_passed_sequences_and_count(self, random_only=False):
        if random_only:
            return {self.cut_adapters(key): value for key, value in self.get_sequences_and_count()}
        else:
            c = Counter(self.get_filter_passed_sequences())
            return c.most_common()

    def filter_function(self, read):
        has_forward = read[: self.forward_adapter_length] == self.forward_adapter \
            or self.forward_adapter_length == 0
        has_reverse = read[-self.reverse_adapter_length:] == self.reverse_adapter \
            or self.reverse_adapter_length == 0
        match_random_region_len = abs(
            len(read) - self.target_length) <= self.tolerance
        return has_forward and has_reverse and match_random_region_len

    def get_filter_passed_sequences(self, random_only=False):
        self.filter_passed = list(filter(self.filter_function, self.raw_reads))
        if random_only:
            return [self.cut_adapters(read) for read in self.filter_passed]
        return self.filter_passed

    def cut_adapters(self, seq):
        if self.reverse_adapter_length == 0:
            ret = seq[self.forward_adapter_length:]
        else:
            ret = seq[self.forward_adapter_length: -
                      self.reverse_adapter_length]
        if self.max_len is not None:
            return ret[len(ret) // 2 - self.max_len // 2: len(ret) // 2 - self.max_len // 2 + self.max_len]
        else:
            return ret

    def __str__(self):
        return f"experiment of {len(self.raw_reads)} raw reads"

    def get_dataloader(self, min_count=1, test_size=0.1, batch_size=512, shuffle=True, use_cuda=True):
        from sklearn.model_selection import train_test_split
        from torch.utils.data import DataLoader

        self.min_count = min_count
        kwargs = {'num_workers': 1, 'pin_memory': True} if (
            use_cuda and torch.cuda.is_available()) else {}
        # load RAPT1-4R and filter reads to count>1, then make it to one hot encoded tensor
        c = self.get_filter_passed_sequences(random_only=True)
        sequences = list(
            filter(lambda seq_count: seq_count[1] >= min_count, Counter(c).most_common()))
        seq, _ = zip(*sequences)

        train_test = np.array(list(map(one_hot_index, seq)))
        logger.info(f"# of sequences -> {len(train_test)}")
        train_data, test_data = train_test_split(
            train_test, test_size=test_size, shuffle=shuffle)
        train_data = torch.from_numpy(train_data).long()
        test_data = torch.from_numpy(test_data).long()
        train_loader = DataLoader(
            train_data, batch_size=batch_size, shuffle=True,  **kwargs)
        test_loader = DataLoader(
            test_data,  batch_size=batch_size, shuffle=False, **kwargs)
        return train_loader, test_loader
```

