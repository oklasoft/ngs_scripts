require 'uri'
module URI
  class O3 < Generic
    if RUBY_VERSION < "2.0"
      USE_REGISTRY = true
    end
    def initialize(*args)
      super(*args)
      if RUBY_VERSION >= "2.0"
        @account = @host
      end
      @host = "o3.omrf.org"
      @port = 443
    end
    if RUBY_VERSION < "2.0"
      def account
        self.registry
      end
    else
      def account
        @account
      end
    end
    def container
      self.path.scan(/^\/([^\/]+)\/(.*)/)[0].first
    end
    def object
      self.path.scan(/^\/([^\/]+)\/(.*)/)[0].last
    end
    def request_uri
      "/v1/#{account()}#{self.path}"
    end
    def os_storage_url
      "https://#{self.host}/v1/#{account()}"
    end
    def https_uri
      URI::HTTPS.build(:host=>self.host,:path=>"/v1/#{account()}#{self.path}")
    end
  end
  @@schemes['O3'] = O3
end
