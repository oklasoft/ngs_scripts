module URI
  class O3 < Generic
    USE_REGISTRY = true
    def initialize(*args)
      super(*args)
      @host = "o3.omrf.org"
      @port = 443
    end
    def account
      self.registry
    end
    def container
      self.path.scan(/^\/([^\/]+)\/(.*)/)[0].first
    end
    def object
      self.path.scan(/^\/([^\/]+)\/(.*)/)[0].last
    end
    def request_uri
      "/v1/#{self.registry}#{self.path}"
    end
    def os_storage_url
      "https://#{self.host}/v1/#{account()}"
    end
    def https_uri
      URI::HTTPS.build(:host=>self.host,:path=>"/v1/#{self.registry}#{self.path}")
    end
  end
  @@schemes['O3'] = O3
end
